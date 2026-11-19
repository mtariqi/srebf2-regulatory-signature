#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase15/phase15c_prepare_inputs.py
# Purpose     : Prepare frozen SREBF2 functional-analysis inputs
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Prepare frozen SREBF2 functional-analysis inputs.

Uses verified primary DESeq2 Wald statistics.
Does not modify upstream analyses.
Refuses to overwrite existing outputs.
"""

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path("results/phase15")
DE_FILE = Path(
    "results/de/discovery_v1/all12/"
    "KO_vs_WT_depletion.tsv"
)
MATRIX_FILE = Path(
    "results/phase12/"
    "SREBF2_HELA_REGULATORY_SIGNAL.tsv"
)

OUTPUTS = {
    "universe": ROOT / "protein_coding_universe.tsv",
    "increased": ROOT / "increased_in_KO.tsv",
    "decreased": ROOT / "decreased_in_KO.tsv",
    "ranked": ROOT / "gsea_primary_wald_ranked.tsv",
    "canonical": ROOT / "canonical_srebf2_genes.tsv",
    "manifest": ROOT / "phase15c_manifest.json",
}

for path in OUTPUTS.values():
    if path.exists():
        raise FileExistsError(
            f"Refusing overwrite: {path}"
        )

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024), b""
        ):
            h.update(chunk)
    return h.hexdigest()

de = pd.read_csv(DE_FILE, sep="\t")
matrix = pd.read_csv(MATRIX_FILE, sep="\t")

assert len(de) == 17158
assert len(matrix) == 17158
assert de["gene_id"].is_unique
assert matrix["gene_id"].is_unique

required = [
    "gene_id", "log2FoldChange", "stat",
    "pvalue", "padj"
]

if not set(required).issubset(de.columns):
    raise ValueError("Required DESeq2 fields missing")

merged = matrix.merge(
    de[
        [
            "gene_id", "log2FoldChange",
            "lfcSE", "stat", "pvalue", "padj"
        ]
    ],
    on="gene_id",
    how="inner",
    validate="one_to_one"
)

assert len(merged) == 17158

delta = np.abs(
    merged["log2FC_primary"].to_numpy()
    - merged["log2FoldChange"].to_numpy()
)

if not np.isfinite(delta).all():
    raise ValueError("Nonfinite log2FC values")

if delta.max() > 1e-8:
    raise ValueError("Primary fold changes disagree")

pc = merged[
    merged["gene_type"] == "protein_coding"
].copy()

assert len(pc) == 12485

if not np.isfinite(
    pc["stat"].to_numpy()
).all():
    raise ValueError("Nonfinite Wald statistics")

increased = pc[
    (pc["signature_selected"] == True)
    & (pc["direction"] == "increased_in_KO")
].copy()

decreased = pc[
    (pc["signature_selected"] == True)
    & (pc["direction"] == "decreased_in_KO")
].copy()

assert len(increased) == 330
assert len(decreased) == 625

columns = [
    "gene_id", "gene_symbol", "baseMean",
    "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"
]

ranked = (
    pc[["gene_id", "gene_symbol", "stat"]]
    .sort_values(
        ["stat", "gene_id"],
        ascending=[False, True]
    )
)

targets = [
    "HMGCR", "HMGCS1", "SQLE",
    "LDLR", "INSIG1", "MSMO1", "CYP51A1"
]

canonical = pc[
    pc["gene_symbol"].isin(targets)
].copy()

assert len(canonical) == 7

ROOT.mkdir(parents=True, exist_ok=True)

tables = {
    "universe": pc[columns].sort_values("gene_id"),
    "increased": increased[columns].sort_values("gene_id"),
    "decreased": decreased[columns].sort_values("gene_id"),
    "ranked": ranked,
    "canonical": canonical[
        columns + [
            "direction", "log2FC_fbs",
            "padj_fbs",
            "log2FC_sensitivity",
            "padj_sensitivity"
        ]
    ].sort_values("gene_symbol"),
}

for name, table in tables.items():
    table.to_csv(
        OUTPUTS[name],
        sep="\t",
        index=False,
        mode="x",
        lineterminator="\n"
    )

manifest = {
    "phase": "15C",
    "primary_contrast": "KO_vs_WT_depletion",
    "ranking": "original DESeq2 signed Wald statistic",
    "gene_identifier": "versioned Ensembl gene ID",
    "protein_coding_universe": len(pc),
    "increased_genes": len(increased),
    "decreased_genes": len(decreased),
    "canonical_genes": len(canonical),
    "de_source_sha256": sha256(DE_FILE),
    "matrix_source_sha256": sha256(MATRIX_FILE),
    "maximum_log2fc_difference": float(delta.max()),
}

with OUTPUTS["manifest"].open("x") as f:
    json.dump(manifest, f, indent=2)
    f.write("\n")

print("=== Phase 15C input preparation ===")

for name, path in OUTPUTS.items():
    print(f"{name}: {path}")
    print("SHA-256:", sha256(path))

print("\nUniverse:", len(pc))
print("Increased:", len(increased))
print("Decreased:", len(decreased))
print("Ranked:", len(ranked))
print("Canonical:", len(canonical))

print("\nPASS: Phase 15C completed")
print("Frozen upstream files unchanged.")
