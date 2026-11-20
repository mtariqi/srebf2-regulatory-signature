#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path("results")
OUT = ROOT / "phase18"
OUT.mkdir(parents=True, exist_ok=True)

CORE = [
    "TM7SF2", "CYP51A1", "MSMO1", "SQLE",
    "FDFT1", "HSD17B7", "LSS", "NSDHL",
    "DHCR7", "EBP", "SC5D", "DHCR24"
]

def read(path):
    return pd.read_csv(path, sep="\t")

# --------------------------------------------------
# 1. Frozen discovery reference
# --------------------------------------------------

annotations = read(
    ROOT / "phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
)

primary = read(
    ROOT / "de/discovery_v1/all12/KO_vs_WT_depletion.tsv"
)

core = annotations.loc[
    annotations["gene_symbol"].isin(CORE),
    ["gene_id", "gene_symbol"]
].copy()

assert len(core) == 12
assert core["gene_id"].is_unique
assert core["gene_symbol"].is_unique

core = core.merge(
    primary[["gene_id", "log2FoldChange", "padj"]],
    on="gene_id",
    validate="one_to_one"
)

core = core.rename(columns={
    "log2FoldChange": "discovery_log2FC",
    "padj": "discovery_padj"
})

assert (core["discovery_log2FC"] < 0).all()
assert core["discovery_log2FC"].notna().all()

core["unversioned_id"] = (
    core["gene_id"].str.replace(
        r"\.\d+$", "", regex=True
    )
)

assert core["unversioned_id"].is_unique

# --------------------------------------------------
# 2. HAP1: two separate serum conditions
# --------------------------------------------------

records = []

hap1_files = {
    "HAP1_minus_FBS":
        ROOT / "validation/GSE147770/interaction_v1/"
        "ko_vs_wt_minus_FBS_v1.tsv",
    "HAP1_plus_FBS":
        ROOT / "validation/GSE147770/interaction_v1/"
        "ko_vs_wt_plus_FBS_v1.tsv"
}

for condition, filename in hap1_files.items():

    df = read(filename)

    assert df["unversioned_id"].is_unique
    assert df["contrast"].nunique() == 1

    merged = core.merge(
        df[[
            "unversioned_id",
            "log2FoldChange",
            "padj"
        ]],
        on="unversioned_id",
        how="left",
        validate="one_to_one"
    )

    for row in merged.itertuples(index=False):
        records.append({
            "dataset": "GSE147770",
            "comparison": condition,
            "gene_id": row.gene_id,
            "gene_symbol": row.gene_symbol,
            "discovery_log2FC": row.discovery_log2FC,
            "discovery_padj": row.discovery_padj,
            "validation_log2FC": row.log2FoldChange,
            "validation_padj": row.padj,
            "orientation_status": "KO_vs_WT_labelled"
        })

# --------------------------------------------------
# 3. HepG2: three separate shRNA contrasts
# --------------------------------------------------

mapping = read(
    ROOT / "validation/GSE271000/processed_DE_v1/"
    "VALIDATION_GENE_MAPPING.tsv"
)

manifest_path = (
    ROOT / "validation/GSE271000/processed_DE_v1/"
    "VALIDATION_MANIFEST.json"
)

manifest = json.loads(manifest_path.read_text())

orientation_verified = (
    manifest.get("contrast_orientation_verified") is True
)

mapping["mapping_eligible"] = (
    mapping["mapping_eligible"]
    .astype(str).str.lower().eq("true")
)

mapping["hairpin"] = (
    mapping["hairpin"].astype(str)
)

for hairpin in ["664", "665", "667"]:

    sub = mapping.loc[
        mapping["hairpin"] == hairpin
    ].copy()

    if sub["gene_id"].duplicated().any():
        raise ValueError(
            f"Duplicate discovery IDs for hairpin {hairpin}"
        )

    sub.loc[
        ~sub["mapping_eligible"],
        ["log2FoldChange", "validation_padj"]
    ] = np.nan

    merged = core.merge(
        sub[[
            "gene_id",
            "gene_symbol",
            "log2FoldChange",
            "validation_padj"
        ]],
        on="gene_id",
        how="left",
        suffixes=("", "_validation"),
        validate="one_to_one"
    )

    mismatched = (
        merged["gene_symbol_validation"].notna() &
        (
            merged["gene_symbol"] !=
            merged["gene_symbol_validation"]
        )
    )

    if mismatched.any():
        raise ValueError(
            f"Symbol mismatch for hairpin {hairpin}"
        )

    for row in merged.itertuples(index=False):
        records.append({
            "dataset": "GSE271000",
            "comparison": f"HepG2_shRNA_{hairpin}",
            "gene_id": row.gene_id,
            "gene_symbol": row.gene_symbol,
            "discovery_log2FC": row.discovery_log2FC,
            "discovery_padj": row.discovery_padj,
            "validation_log2FC": row.log2FoldChange,
            "validation_padj": row.validation_padj,
            "orientation_status": (
                "verified" if orientation_verified
                else "UNVERIFIED"
            )
        })

# --------------------------------------------------
# 4. Per-gene concordance
# --------------------------------------------------

out = pd.DataFrame(records)

out["validation_log2FC"] = pd.to_numeric(
    out["validation_log2FC"],
    errors="coerce"
)

out["validation_padj"] = pd.to_numeric(
    out["validation_padj"],
    errors="coerce"
)

out["evaluable"] = (
    np.isfinite(out["validation_log2FC"]) &
    out["validation_log2FC"].ne(0)
)

out["negative_validation_LFC"] = (
    out["evaluable"] &
    out["validation_log2FC"].lt(0)
)

out["direction_concordant"] = (
    out["evaluable"] &
    (
        np.sign(out["discovery_log2FC"]) ==
        np.sign(out["validation_log2FC"])
    )
)

out["negative_and_FDR05"] = (
    out["negative_validation_LFC"] &
    out["validation_padj"].lt(0.05)
)

out["direction_status"] = np.select(
    [
        ~out["evaluable"],
        out["direction_concordant"]
    ],
    [
        "NOT_EVALUABLE",
        "CONCORDANT"
    ],
    default="DISCORDANT"
)

assert len(out) == 60

out.to_csv(
    OUT / "phase18c_sterol_core_gene_level.tsv",
    sep="\t",
    index=False
)

# --------------------------------------------------
# 5. Per-contrast summaries
# --------------------------------------------------

summary = []

for comparison, group in out.groupby(
    "comparison", sort=False
):

    evaluable = group.loc[group["evaluable"]]
    n = len(evaluable)

    concordant = int(
        evaluable["direction_concordant"].sum()
    )

    if n >= 3:
        rho, rho_p = spearmanr(
            evaluable["discovery_log2FC"],
            evaluable["validation_log2FC"]
        )
    else:
        rho, rho_p = np.nan, np.nan

    summary.append({
        "comparison": comparison,
        "dataset": group["dataset"].iloc[0],
        "orientation_status":
            group["orientation_status"].iloc[0],
        "core_genes": 12,
        "evaluable_genes": n,
        "direction_concordant": concordant,
        "direction_concordance_fraction":
            concordant / n if n else np.nan,
        "negative_FDR05": int(
            group["negative_and_FDR05"].sum()
        ),
        "spearman_rho": rho,
        "spearman_pvalue_exploratory": rho_p,
        "binomial_pvalue_exploratory": (
            binomtest(
                concordant, n, 0.5
            ).pvalue if n else np.nan
        )
    })

summary = pd.DataFrame(summary)

summary.to_csv(
    OUT / "phase18c_sterol_core_summary.tsv",
    sep="\t",
    index=False
)

print("=== PHASE 18C: STEROL-CORE VALIDATION ===")
print(summary.to_string(index=False))

print("\n=== PER-GENE VALIDATION FOLD CHANGES ===")
print(
    out.pivot(
        index="gene_symbol",
        columns="comparison",
        values="validation_log2FC"
    ).round(3).to_string()
)

print("\n=== IMPORTANT LIMITATION ===")
print(
    "HepG2 orientation verified:",
    orientation_verified
)
print(
    "HepG2 directional findings remain provisional "
    "until the deposited contrast definition is confirmed."
)

print("\nPASS: Phase 18C analysis completed")
