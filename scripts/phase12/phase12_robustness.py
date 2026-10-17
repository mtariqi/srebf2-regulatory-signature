#!/usr/bin/env python3
"""Phase 12J: Robustness of direction-specific HeLa regulatory associations."""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from pathlib import Path

INPUT = Path("results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv")
OUTPUT = Path("results/phase12/PHASE12_ROBUSTNESS_MODELS.tsv")
VIF_OUTPUT = Path("results/phase12/PHASE12_VIF.tsv")

df = pd.read_csv(INPUT, sep="\t", low_memory=False)
df = df.loc[df["_merge"] == "both"].copy()

df["signature"] = (
    df["signature_selected"].astype(str).str.lower() == "true"
)

df["log_baseMean"] = np.log10(df["baseMean"] + 1)
df["log_gene_length"] = np.log10(df["gene_length"])

df["biotype_group"] = df["gene_type"].where(
    df["gene_type"].isin(["protein_coding", "lncRNA"]),
    "other"
)

tracks = {
    "H3K27ac": "mean_promoter_H3K27ac",
    "SREBP2_CholDep": "mean_promoter_SREBP2_CholDep",
    "Input_CholDep": "mean_promoter_Input_CholDep",
}

for name, column in tracks.items():
    values = np.log1p(
        pd.to_numeric(df[column], errors="coerce")
    )
    df[f"z_{name}"] = (
        values - values.mean()
    ) / values.std()

required = [
    "log_baseMean",
    "log_gene_length",
    "biotype_group",
    "z_H3K27ac",
    "z_SREBP2_CholDep",
    "z_Input_CholDep",
]

df = df.dropna(subset=required).copy()

print("========== PHASE 12J ==========")
print("Matched genes:", len(df))
print("Signature genes:", int(df["signature"].sum()))

# --------------------------------------------------
# VIF analysis
# --------------------------------------------------

print("\n========== VIF DIAGNOSTICS ==========")

predictors = [
    "z_H3K27ac",
    "z_SREBP2_CholDep",
    "z_Input_CholDep",
    "log_baseMean",
    "log_gene_length",
]

X = df[predictors].copy()
X = (X - X.mean()) / X.std()
X.insert(0, "intercept", 1.0)

vif_results = []

for index, predictor in enumerate(X.columns):
    if predictor == "intercept":
        continue

    vif_results.append({
        "predictor": predictor,
        "VIF": variance_inflation_factor(
            X.values, index
        ),
    })

vif_df = pd.DataFrame(vif_results)
vif_df.to_csv(VIF_OUTPUT, sep="\t", index=False)

print(vif_df.to_string(index=False))

# --------------------------------------------------
# Define direction-specific comparisons
# --------------------------------------------------

comparisons = {
    "decreased_vs_background": df.loc[
        (~df["signature"]) |
        (df["direction"] == "decreased_in_KO")
    ].copy(),

    "increased_vs_background": df.loc[
        (~df["signature"]) |
        (df["direction"] == "increased_in_KO")
    ].copy(),

    "increased_vs_decreased": df.loc[
        df["signature"]
    ].copy(),
}

specifications = {
    "base": (
        "outcome ~ z_H3K27ac + z_SREBP2_CholDep"
        " + log_baseMean + log_gene_length"
        " + C(biotype_group)"
    ),

    "input_adjusted": (
        "outcome ~ z_H3K27ac + z_SREBP2_CholDep"
        " + z_Input_CholDep"
        " + log_baseMean + log_gene_length"
        " + C(biotype_group)"
    ),
}

results = []

for comparison, subset in comparisons.items():

    if comparison == "increased_vs_decreased":
        subset["outcome"] = (
            subset["direction"] == "increased_in_KO"
        ).astype(int)
    else:
        subset["outcome"] = subset["signature"].astype(int)

    for population, data in [
        ("all_biotypes", subset),
        (
            "protein_coding",
            subset.loc[
                subset["gene_type"] == "protein_coding"
            ].copy()
        ),
    ]:

        for specification, formula in specifications.items():

            if population == "protein_coding":
                formula = formula.replace(
                    " + C(biotype_group)", ""
                )

            fit = smf.logit(
                formula,
                data=data
            ).fit(
                disp=False,
                maxiter=200
            )

            confidence = fit.conf_int()

            for predictor in [
                "z_H3K27ac",
                "z_SREBP2_CholDep",
            ]:

                results.append({
                    "comparison": comparison,
                    "population": population,
                    "specification": specification,
                    "predictor": predictor,
                    "n": int(fit.nobs),
                    "cases": int(data["outcome"].sum()),
                    "odds_ratio": np.exp(
                        fit.params[predictor]
                    ),
                    "ci_lower": np.exp(
                        confidence.loc[predictor, 0]
                    ),
                    "ci_upper": np.exp(
                        confidence.loc[predictor, 1]
                    ),
                    "p_value": fit.pvalues[predictor],
                    "converged": fit.mle_retvals["converged"],
                })

            print(
                f"{comparison} | {population} | "
                f"{specification} | "
                f"n={int(fit.nobs)} | "
                f"converged={fit.mle_retvals['converged']}"
            )

result = pd.DataFrame(results)
result.to_csv(OUTPUT, sep="\t", index=False)

print("\n========== ROBUSTNESS RESULTS ==========")
print(result.to_string(index=False))

print("\nSaved:", OUTPUT)
print("Saved:", VIF_OUTPUT)
print("\nPHASE 12J COMPLETE")
