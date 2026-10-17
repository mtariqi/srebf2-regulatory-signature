#!/usr/bin/env python3
"""Phase 12I: Direction-specific adjusted regulatory-signal models."""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

INPUT = "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
OUTPUT = "results/phase12/PHASE12_DIRECTION_MODELS.tsv"

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

signals = {
    "H3K27ac": "mean_promoter_H3K27ac",
    "SREBP2_FBS": "mean_promoter_SREBP2_FBS",
    "SREBP2_CholDep": "mean_promoter_SREBP2_CholDep",
}

for name, column in signals.items():
    values = np.log1p(pd.to_numeric(df[column], errors="coerce"))
    df[f"z_{name}"] = (values - values.mean()) / values.std()

required = [
    "log_baseMean", "log_gene_length", "biotype_group",
    "z_H3K27ac", "z_SREBP2_FBS", "z_SREBP2_CholDep"
]

df = df.dropna(subset=required).copy()

covariates = (
    "log_baseMean + log_gene_length + C(biotype_group)"
)

model_specs = {
    "H3K27ac": f"outcome ~ z_H3K27ac + {covariates}",
    "SREBP2_FBS": f"outcome ~ z_SREBP2_FBS + {covariates}",
    "SREBP2_CholDep": f"outcome ~ z_SREBP2_CholDep + {covariates}",
    "Combined_CholDep":
        f"outcome ~ z_H3K27ac + z_SREBP2_CholDep + {covariates}",
}

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

results = []

for comparison, subset in comparisons.items():

    if comparison == "increased_vs_decreased":
        subset["outcome"] = (
            subset["direction"] == "increased_in_KO"
        ).astype(int)
    else:
        subset["outcome"] = subset["signature"].astype(int)

    print(
        f"\n{comparison}: "
        f"n={len(subset)}, cases={subset['outcome'].sum()}"
    )

    for model_name, formula in model_specs.items():
        fit = smf.logit(formula, data=subset).fit(
            disp=False, maxiter=200
        )

        ci = fit.conf_int()

        for predictor in fit.params.index:
            if not predictor.startswith("z_"):
                continue

            results.append({
                "comparison": comparison,
                "model": model_name,
                "predictor": predictor,
                "n": int(fit.nobs),
                "odds_ratio": np.exp(fit.params[predictor]),
                "ci_lower": np.exp(ci.loc[predictor, 0]),
                "ci_upper": np.exp(ci.loc[predictor, 1]),
                "p_value": fit.pvalues[predictor],
                "converged": fit.mle_retvals["converged"],
            })

result = pd.DataFrame(results)

result["fdr_bh"] = multipletests(
    result["p_value"], method="fdr_bh"
)[1]

result.to_csv(OUTPUT, sep="\t", index=False)

print("\n========== PHASE 12I RESULTS ==========")
print(result.to_string(index=False))
print("\nSaved:", OUTPUT)
