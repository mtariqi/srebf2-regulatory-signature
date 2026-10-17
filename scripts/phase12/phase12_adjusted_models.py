#!/usr/bin/env python3
"""Phase 12H: Adjusted HeLa regulatory-signal associations."""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

INPUT = "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
OUTPUT = "results/phase12/PHASE12_ADJUSTED_MODELS.tsv"

df = pd.read_csv(INPUT, sep="\t", low_memory=False)
df = df.loc[df["_merge"] == "both"].copy()

df["signature"] = (
    df["signature_selected"].astype(str).str.lower() == "true"
).astype(int)

df["log_baseMean"] = np.log10(df["baseMean"] + 1)
df["log_gene_length"] = np.log10(df["gene_length"])

# Group sparse biotypes.
df["biotype_group"] = df["gene_type"].where(
    df["gene_type"].isin(["protein_coding", "lncRNA"]),
    "other"
)

signals = {
    "H3K27ac": "mean_promoter_H3K27ac",
    "SREBP2_FBS": "mean_promoter_SREBP2_FBS",
    "SREBP2_CholDep": "mean_promoter_SREBP2_CholDep",
}

for label, column in signals.items():
    transformed = np.log1p(
        pd.to_numeric(df[column], errors="coerce")
    )
    df[f"z_{label}"] = (
        (transformed - transformed.mean()) / transformed.std()
    )

# Use the same complete-case population in every model.
required = [
    "signature", "log_baseMean", "log_gene_length",
    "biotype_group",
    "z_H3K27ac", "z_SREBP2_FBS", "z_SREBP2_CholDep",
]

df = df.dropna(subset=required).copy()

covariates = (
    "log_baseMean + log_gene_length + C(biotype_group)"
)

models = {
    "M1_H3K27ac":
        f"signature ~ z_H3K27ac + {covariates}",
    "M2_SREBP2_FBS":
        f"signature ~ z_SREBP2_FBS + {covariates}",
    "M3_SREBP2_CholDep":
        f"signature ~ z_SREBP2_CholDep + {covariates}",
    "M4_H3K27ac_FBS":
        f"signature ~ z_H3K27ac + z_SREBP2_FBS + {covariates}",
    "M5_H3K27ac_CholDep":
        f"signature ~ z_H3K27ac + z_SREBP2_CholDep + {covariates}",
}

results = []

def fit_model(name, formula, data):
    fit = smf.logit(formula, data=data).fit(
        disp=False, maxiter=200
    )

    confidence = fit.conf_int()

    for variable in fit.params.index:
        if not variable.startswith("z_"):
            continue

        results.append({
            "model": name,
            "predictor": variable,
            "n": int(fit.nobs),
            "odds_ratio": np.exp(fit.params[variable]),
            "ci_lower": np.exp(confidence.loc[variable, 0]),
            "ci_upper": np.exp(confidence.loc[variable, 1]),
            "p_value": fit.pvalues[variable],
            "aic": fit.aic,
        })

    print(f"{name}: converged={fit.mle_retvals['converged']}")
    return fit


print("========== PRIMARY MODELS ==========")
print("Analysis genes:", len(df))
print("Signature genes:", df["signature"].sum())

for name, formula in models.items():
    fit_model(name, formula, df)

print("\n========== DIRECTION-SPECIFIC MODELS ==========")

print("Observed direction labels:")
print(df.loc[df["signature"] == 1, "direction"].value_counts())

# Inspect the labels before choosing direction-specific subsets.
direction_values = df.loc[
    df["signature"] == 1, "direction"
].dropna().unique()

print("Unique signature directions:", list(direction_values))

result = pd.DataFrame(results)

result["fdr_bh"] = multipletests(
    result["p_value"],
    method="fdr_bh"
)[1]

result.to_csv(OUTPUT, sep="\t", index=False)

print("\n========== ADJUSTED RESULTS ==========")
print(result.to_string(index=False))
print("\nSaved:", OUTPUT)
