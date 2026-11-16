#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase13/phase13e_adjusted_enrichment.py
# Purpose     : Phase 13E: exploratory covariate-adjusted E117 state associations
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 13E: exploratory covariate-adjusted E117 state associations."""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

INPUT = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)

OUTPUT = Path(
    "results/phase13/SREBF2_E117_ADJUSTED_ENRICHMENT.tsv"
)

STATES = [
    "1_TssA",
    "2_TssAFlnk",
    "5_TxWk",
    "7_Enh",
    "14_ReprPCWk",
    "15_Quies",
]

CONTRASTS = [
    ("decreased_in_KO", "background"),
    ("increased_in_KO", "background"),
    ("increased_in_KO", "decreased_in_KO"),
]

EXPECTED = {
    "background": 15781,
    "decreased_in_KO": 789,
    "increased_in_KO": 527,
}

if OUTPUT.exists():
    raise FileExistsError(f"Refusing to overwrite {OUTPUT}")

df = pd.read_csv(INPUT, sep="\t")

if df["gene_id"].duplicated().any():
    raise ValueError("Duplicate gene IDs")

actual = df["signature_group"].value_counts().to_dict()
if actual != EXPECTED:
    raise ValueError(f"Unexpected signature counts: {actual}")

# Unannotated TSSs are missing reference-state observations.
df = df[df["chromhmm_state"] != "UNANNOTATED"].copy()

df["baseMean"] = pd.to_numeric(df["baseMean"], errors="coerce")
df["gene_length"] = pd.to_numeric(df["gene_length"], errors="coerce")

if (df["baseMean"].dropna() < 0).any():
    raise ValueError("Negative baseline expression detected")

if (df["gene_length"].dropna() <= 0).any():
    raise ValueError("Nonpositive gene length detected")

df["log_baseMean"] = np.log10(1 + df["baseMean"])
df["log_gene_length"] = np.log10(df["gene_length"])

def biotype_group(value):
    if value == "protein_coding":
        return "protein_coding"
    if value == "lncRNA":
        return "lncRNA"
    return "other"

df["biotype_group"] = df["gene_type"].map(biotype_group)

required = [
    "log_baseMean",
    "log_gene_length",
    "biotype_group",
]

results = []

for group_a, group_b in CONTRASTS:

    sub = df[
        df["signature_group"].isin([group_a, group_b])
    ].copy()

    n_before = len(sub)
    sub = sub.dropna(subset=required)

    sub["outcome"] = (
        sub["signature_group"] == group_a
    ).astype(int)

    if sub["outcome"].nunique() != 2:
        raise ValueError(f"Invalid contrast: {group_a}, {group_b}")

    for state in STATES:

        data = sub.copy()

        data["state_present"] = (
            data["chromhmm_state"] == state
        ).astype(int)

        if data["state_present"].nunique() != 2:
            raise ValueError(
                f"No variation for {state} in {group_a} vs {group_b}"
            )

        formula = (
            "outcome ~ state_present + log_baseMean "
            "+ log_gene_length + C(biotype_group)"
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")

            model = smf.glm(
                formula=formula,
                data=data,
                family=sm.families.Binomial(),
            ).fit(maxiter=200)

        if not model.converged:
            raise RuntimeError(
                f"Nonconverged model: {group_a} vs {group_b}, {state}"
            )

        beta = float(model.params["state_present"])
        se = float(model.bse["state_present"])
        pvalue = float(model.pvalues["state_present"])

        if not np.all(np.isfinite([beta, se, pvalue])):
            raise RuntimeError(
                f"Nonfinite estimate: {group_a} vs {group_b}, {state}"
            )

        ci = model.conf_int().loc["state_present"]

        results.append({
            "contrast": f"{group_a}_vs_{group_b}",
            "state": state,
            "n_total": len(data),
            "n_dropped_covariates": n_before - len(data),
            "n_group_a": int(data["outcome"].sum()),
            "n_group_b": int((1 - data["outcome"]).sum()),
            "n_state_group_a": int(
                ((data["outcome"] == 1) &
                 (data["state_present"] == 1)).sum()
            ),
            "n_state_group_b": int(
                ((data["outcome"] == 0) &
                 (data["state_present"] == 1)).sum()
            ),
            "adjusted_or": float(np.exp(beta)),
            "ci95_lower": float(np.exp(ci.iloc[0])),
            "ci95_upper": float(np.exp(ci.iloc[1])),
            "pvalue": pvalue,
            "converged": bool(model.converged),
            "warnings": "; ".join(
                sorted(set(str(w.message) for w in caught))
            ),
        })

if len(results) != 18:
    raise ValueError(f"Expected 18 models, obtained {len(results)}")

qvalues = multipletests(
    [r["pvalue"] for r in results],
    method="fdr_bh",
)[1]

for r, q in zip(results, qvalues):
    r["fdr_bh_18tests"] = float(q)

out = pd.DataFrame(results)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPUT, sep="\t", index=False, mode="x")

print("=== Phase 13E adjusted ChromHMM enrichment ===")
print("Status: PASS")
print("Models:", len(out))
print("Converged:", int(out["converged"].sum()))
print("FDR < 0.05:", int((out["fdr_bh_18tests"] < 0.05).sum()))

print("\n=== FDR-significant associations ===")

sig = out[out["fdr_bh_18tests"] < 0.05].sort_values(
    "fdr_bh_18tests"
)

for _, row in sig.iterrows():
    print(
        f'{row["contrast"]} | {row["state"]} | '
        f'OR={row["adjusted_or"]:.3f} | '
        f'95% CI={row["ci95_lower"]:.3f}-'
        f'{row["ci95_upper"]:.3f} | '
        f'FDR={row["fdr_bh_18tests"]:.3e}'
    )

print("\nOutput:", OUTPUT)
