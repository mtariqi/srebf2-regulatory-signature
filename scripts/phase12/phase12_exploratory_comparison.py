#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

INPUT = "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
OUTPUT = "results/phase12/PHASE12_EXPLORATORY_COMPARISONS.tsv"

df = pd.read_csv(INPUT, sep="\t", low_memory=False)

# Explicitly retain only genes with exact annotation matches.
df = df.loc[df["_merge"] == "both"].copy()

if df["signature_selected"].dtype != bool:
    df["signature_selected"] = (
        df["signature_selected"].astype(str).str.lower() == "true"
    )

signals = [
    "mean_promoter_H3K27ac",
    "mean_genebody_H3K27ac",
    "mean_promoter_SREBP2_FBS",
    "mean_promoter_SREBP2_CholDep",
    "mean_promoter_Input_CholDep",
]

results = []

for signal in signals:
    signature = pd.to_numeric(
        df.loc[df["signature_selected"], signal],
        errors="coerce"
    ).dropna()

    background = pd.to_numeric(
        df.loc[~df["signature_selected"], signal],
        errors="coerce"
    ).dropna()

    test = mannwhitneyu(
        signature,
        background,
        alternative="two-sided",
        method="asymptotic"
    )

    n1 = len(signature)
    n0 = len(background)

    # Rank-biserial correlation: positive means higher
    # signal in the signature group.
    rank_biserial = 2 * test.statistic / (n1 * n0) - 1

    results.append({
        "signal": signal,
        "n_signature": n1,
        "n_background": n0,
        "median_signature": signature.median(),
        "median_background": background.median(),
        "median_difference": signature.median() - background.median(),
        "rank_biserial": rank_biserial,
        "p_value": test.pvalue,
    })

result = pd.DataFrame(results)

# Benjamini-Hochberg adjustment.
p = result["p_value"].to_numpy()
order = np.argsort(p)
rank = np.arange(1, len(p) + 1)
adjusted = np.minimum.accumulate(
    (p[order] * len(p) / rank)[::-1]
)[::-1]

q = np.empty_like(adjusted)
q[order] = np.minimum(adjusted, 1.0)

result["fdr_bh"] = q

result.to_csv(OUTPUT, sep="\t", index=False)

print("\nPHASE 12G — EXPLORATORY SIGNAL COMPARISONS")
print("Genes analyzed:", len(df))
print("Signature genes:", int(df["signature_selected"].sum()))
print("Non-signature genes:", int((~df["signature_selected"]).sum()))
print()
print(result.to_string(index=False))
print("\nSaved:", OUTPUT)
