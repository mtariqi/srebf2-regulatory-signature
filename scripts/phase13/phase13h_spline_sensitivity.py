#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase13/phase13h_spline_sensitivity.py
# Purpose     : Phase 13H: flexible covariate adjustment for E117 state enrichment
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 13H: flexible covariate adjustment for E117 state enrichment."""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

SOURCE = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)
OUTPUT = Path(
    "results/phase13/SREBF2_E117_SPLINE_SENSITIVITY.tsv"
)

STATES = [
    "1_TssA", "2_TssAFlnk", "5_TxWk",
    "7_Enh", "14_ReprPCWk", "15_Quies",
]

CONTRASTS = [
    ("decreased_in_KO", "background"),
    ("increased_in_KO", "background"),
    ("increased_in_KO", "decreased_in_KO"),
]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

df = pd.read_csv(SOURCE, sep="\t")

if len(df) != 17097 or df["gene_id"].duplicated().any():
    raise ValueError("Unexpected gene universe or duplicate IDs")

expected = {
    "background": 15781,
    "decreased_in_KO": 789,
    "increased_in_KO": 527,
}

if df["signature_group"].value_counts().to_dict() != expected:
    raise ValueError("Frozen signature groups do not match expected counts")

df = df[df["chromhmm_state"] != "UNANNOTATED"].copy()

for column in ("baseMean", "gene_length"):
    df[column] = pd.to_numeric(df[column], errors="coerce")

if df[["baseMean", "gene_length"]].isna().any().any():
    raise ValueError("Missing covariates")

if not np.isfinite(
    df[["baseMean", "gene_length"]].to_numpy()
).all():
    raise ValueError("Nonfinite covariates")

if (df["baseMean"] < 0).any() or (df["gene_length"] <= 0).any():
    raise ValueError("Invalid expression or gene length")

df["log_baseMean"] = np.log10(df["baseMean"] + 1)
df["log_gene_length"] = np.log10(df["gene_length"])

df["biotype_group"] = np.select(
    [
        df["gene_type"] == "protein_coding",
        df["gene_type"] == "lncRNA",
    ],
    ["protein_coding", "lncRNA"],
    default="other",
)

results = []

for group_a, group_b in CONTRASTS:
    sub = df[
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

        formula = (
            "outcome ~ state_present "
            "+ cr(log_baseMean, df=4) "
            "+ cr(log_gene_length, df=4) "
            "+ C(biotype_group)"
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = smf.glm(
                formula,
                data=data,
                family=sm.families.Binomial(),
            ).fit(maxiter=200)

        if not fit.converged:
            raise RuntimeError(
                f"Nonconverged: {group_a} vs {group_b}, {state}"
            )

        beta = float(fit.params["state_present"])
        p = float(fit.pvalues["state_present"])
        ci = fit.conf_int().loc["state_present"]

        if not np.isfinite(
            [beta, p, ci.iloc[0], ci.iloc[1]]
        ).all():
            raise RuntimeError("Nonfinite model estimate")

        results.append({
            "contrast": f"{group_a}_vs_{group_b}",
            "state": state,
            "n_total": len(data),
            "adjusted_or": float(np.exp(beta)),
            "ci95_lower": float(np.exp(ci.iloc[0])),
            "ci95_upper": float(np.exp(ci.iloc[1])),
            "pvalue": p,
            "converged": fit.converged,
            "warnings": "; ".join(
                sorted(set(str(w.message) for w in caught))
            ),
        })

if len(results) != 18:
    raise RuntimeError("Expected 18 models")

qvalues = multipletests(
    [r["pvalue"] for r in results],
    method="fdr_bh"
)[1]

for r, q in zip(results, qvalues):
    r["fdr_bh_18tests"] = float(q)

out = pd.DataFrame(results)
out.to_csv(OUTPUT, sep="\t", index=False, mode="x")

print("=== Phase 13H spline sensitivity ===")
print("Models:", len(out))
print("Converged:", int(out["converged"].sum()))
print("FDR < 0.05:", int((out["fdr_bh_18tests"] < 0.05).sum()))
print("Models with warnings:",
      int(out["warnings"].fillna("").str.len().gt(0).sum()))

print("\n=== Two cross-model-consistent associations ===")

targets = [
    ("increased_in_KO_vs_background", "15_Quies"),
    ("decreased_in_KO_vs_background", "14_ReprPCWk"),
]

for contrast, state in targets:
    row = out[
        (out["contrast"] == contrast) &
        (out["state"] == state)
    ].iloc[0]

    print(
        f"{contrast} | {state} | "
        f"OR={row.adjusted_or:.3f} | "
        f"95% CI={row.ci95_lower:.3f}-"
        f"{row.ci95_upper:.3f} | "
        f"FDR={row.fdr_bh_18tests:.3e}"
    )

print("\nOutput:", OUTPUT)
