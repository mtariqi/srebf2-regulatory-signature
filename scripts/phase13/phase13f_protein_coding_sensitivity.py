#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase13/phase13f_protein_coding_sensitivity.py
# Purpose     : Phase 13F: protein-coding-only E117 chromatin-state sensitivity
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 13F: protein-coding-only E117 chromatin-state sensitivity."""

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
    "results/phase13/SREBF2_E117_PROTEIN_CODING_SENSITIVITY.tsv"
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
    "background": 11461,
    "decreased_in_KO": 623,
    "increased_in_KO": 328,
}

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

df = pd.read_csv(INPUT, sep="\t")

if df["gene_id"].duplicated().any():
    raise ValueError("Duplicate gene IDs")

df = df.loc[
    (df["gene_type"] == "protein_coding") &
    (df["chromhmm_state"] != "UNANNOTATED")
].copy()

actual = df["signature_group"].value_counts().to_dict()

if actual != EXPECTED:
    raise ValueError(
        f"Unexpected protein-coding group counts: {actual}"
    )

df["baseMean"] = pd.to_numeric(
    df["baseMean"], errors="coerce"
)

df["gene_length"] = pd.to_numeric(
    df["gene_length"], errors="coerce"
)

if df["baseMean"].isna().any() or df["gene_length"].isna().any():
    raise ValueError("Missing numeric covariates")

if not np.isfinite(
    df[["baseMean", "gene_length"]].to_numpy()
).all():
    raise ValueError("Nonfinite numeric covariates")

if (df["baseMean"] < 0).any():
    raise ValueError("Negative baseline expression")

if (df["gene_length"] <= 0).any():
    raise ValueError("Nonpositive gene length")

df["log_baseMean"] = np.log10(1 + df["baseMean"])
df["log_gene_length"] = np.log10(df["gene_length"])

results = []

for group_a, group_b in CONTRASTS:

    sub = df.loc[
        df["signature_group"].isin([group_a, group_b])
    ].copy()

    sub["outcome"] = (
        sub["signature_group"] == group_a
    ).astype(int)

    for state in STATES:

        data = sub.copy()

        data["state_present"] = (
            data["chromhmm_state"] == state
        ).astype(int)

        a = int(
            ((data["outcome"] == 1) &
             (data["state_present"] == 1)).sum()
        )
        c = int(
            ((data["outcome"] == 0) &
             (data["state_present"] == 1)).sum()
        )

        if a == 0 or c == 0:
            raise ValueError(
                f"Sparse zero cell: {group_a} vs {group_b}, {state}"
            )

        formula = (
            "outcome ~ state_present + "
            "log_baseMean + log_gene_length"
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")

            model = smf.glm(
                formula=formula,
                data=data,
                family=sm.families.Binomial()
            ).fit(maxiter=200)

        if not model.converged:
            raise RuntimeError(
                f"Nonconverged model: {group_a}, {group_b}, {state}"
            )

        beta = float(model.params["state_present"])
        se = float(model.bse["state_present"])
        pvalue = float(model.pvalues["state_present"])

        ci = model.conf_int().loc["state_present"]

        if not np.isfinite(
            [beta, se, pvalue, ci.iloc[0], ci.iloc[1]]
        ).all():
            raise RuntimeError(
                f"Nonfinite model estimate: {group_a}, {state}"
            )

        results.append({
            "contrast": f"{group_a}_vs_{group_b}",
            "state": state,
            "n_total": len(data),
            "n_group_a": int(data["outcome"].sum()),
            "n_group_b": int((1-data["outcome"]).sum()),
            "n_state_group_a": a,
            "n_state_group_b": c,
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
    raise ValueError(
        f"Expected 18 models, found {len(results)}"
    )

qvalues = multipletests(
    [r["pvalue"] for r in results],
    method="fdr_bh",
)[1]

for record, qvalue in zip(results, qvalues):
    record["fdr_bh_18tests"] = float(qvalue)

output = pd.DataFrame(results)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
output.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    mode="x"
)

print("=== Phase 13F protein-coding sensitivity ===")
print("Status: PASS")
print("Eligible protein-coding genes:", len(df))
print("Models:", len(output))
print("Converged:", int(output["converged"].sum()))
print(
    "FDR < 0.05:",
    int((output["fdr_bh_18tests"] < 0.05).sum())
)

print("\n=== FDR-significant associations ===")

significant = output.loc[
    output["fdr_bh_18tests"] < 0.05
].sort_values("fdr_bh_18tests")

for _, row in significant.iterrows():
    print(
        f'{row["contrast"]} | {row["state"]} | '
        f'OR={row["adjusted_or"]:.3f} | '
        f'95% CI={row["ci95_lower"]:.3f}-'
        f'{row["ci95_upper"]:.3f} | '
        f'FDR={row["fdr_bh_18tests"]:.3e}'
    )

print("\nModels with warnings:",
      int(output["warnings"].fillna("").str.len().gt(0).sum()))

print("Output:", OUTPUT)
