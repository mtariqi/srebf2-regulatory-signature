#!/usr/bin/env python3

"""Phase 12F: Quality control of HeLa regulatory bigWig signals."""

import pandas as pd
from pathlib import Path

INPUT = Path(
    "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
)

df = pd.read_csv(INPUT, sep="\t", low_memory=False)

signals = [
    "mean_promoter_Input_CholDep",
    "mean_promoter_H3K27ac",
    "mean_genebody_H3K27ac",
    "mean_promoter_SREBP2_FBS",
    "mean_promoter_SREBP2_CholDep",
]

fractions = [
    col for col in df.columns
    if col.startswith("promoter_defined_fraction_")
    or col.startswith("genebody_defined_fraction_")
]

print("========== PHASE 12F: DATASET ==========")

print("Total genes:", len(df))
print("Unique gene IDs:", df["gene_id"].nunique())

print("\nSignature counts:")
print(df["signature_selected"].value_counts(dropna=False))

print("\n========== SIGNAL DISTRIBUTIONS ==========")

for col in signals:
    s = pd.to_numeric(df[col], errors="coerce")

    print(f"\n{col}")
    print("Nonmissing:", s.notna().sum())
    print("Zero:", int((s == 0).sum()))
    print("Negative:", int((s < 0).sum()))
    print("Median:", s.median())
    print("95th percentile:", s.quantile(0.95))
    print("99th percentile:", s.quantile(0.99))
    print("Maximum:", s.max())

print("\n========== COVERAGE FRACTIONS ==========")

for col in fractions:
    s = pd.to_numeric(df[col], errors="coerce")

    print(f"\n{col}")
    print("Median:", s.median())
    print("Minimum:", s.min())
    print("Below 0.90:", int((s < 0.90).sum()))
    print("Below 0.50:", int((s < 0.50).sum()))

print("\n========== SPEARMAN CORRELATIONS ==========")

numeric = df[signals + ["baseMean"]].apply(
    pd.to_numeric, errors="coerce"
)

print(numeric.corr(method="spearman").round(3).to_string())

print("\n========== ANNOTATION MATCHING ==========")

print(df["_merge"].value_counts(dropna=False))

print("\n========== UNMATCHED GENES BY SIGNATURE ==========")

print(
    pd.crosstab(
        df["_merge"],
        df["signature_selected"],
        dropna=False
    ).to_string()
)

print("\n========== QC COMPLETE ==========")
