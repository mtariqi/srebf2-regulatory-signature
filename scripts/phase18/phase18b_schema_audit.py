#!/usr/bin/env python3
from pathlib import Path
import json
import pandas as pd

files = {
    "HeLa discovery": (
        "results/de/discovery_v1/all12/"
        "KO_vs_WT_depletion.tsv"
    ),
    "HeLa gene annotations": (
        "results/phase12/"
        "SREBF2_HELA_REGULATORY_SIGNAL.tsv"
    ),
    "HAP1 minus FBS": (
        "results/validation/GSE147770/"
        "interaction_v1/ko_vs_wt_minus_FBS_v1.tsv"
    ),
    "HAP1 plus FBS": (
        "results/validation/GSE147770/"
        "interaction_v1/ko_vs_wt_plus_FBS_v1.tsv"
    ),
    "HepG2 gene mapping": (
        "results/validation/GSE271000/"
        "processed_DE_v1/VALIDATION_GENE_MAPPING.tsv"
    ),
}

core = {
    "TM7SF2", "CYP51A1", "MSMO1", "SQLE",
    "FDFT1", "HSD17B7", "LSS", "NSDHL",
    "DHCR7", "EBP", "SC5D", "DHCR24"
}

print("=== PHASE 18B: SCHEMA AUDIT ===")

for label, name in files.items():
    path = Path(name)
    print(f"\n=== {label} ===")

    if not path.exists():
        print("MISSING:", path)
        continue

    df = pd.read_csv(path, sep="\t")
    print("Rows:", len(df))
    print("Columns:", df.columns.tolist())

    symbol_cols = [
        c for c in df.columns
        if c.lower() in {
            "gene_symbol", "symbol", "gene_name",
            "external_gene_name"
        }
    ]

    for col in symbol_cols:
        observed = set(df[col].dropna().astype(str))
        print("Core overlap using", col, ":",
              len(core & observed), "/12")

    print("First record:")
    print(df.head(1).to_string(index=False))

manifest = Path(
    "results/validation/GSE271000/"
    "processed_DE_v1/VALIDATION_MANIFEST.json"
)

if manifest.exists():
    obj = json.loads(manifest.read_text())
    print("\n=== HEPG2 CONTRAST PROVENANCE ===")
    print("Assumed orientation:",
          obj.get("assumed_contrast"))
    print("Orientation verified:",
          obj.get("contrast_orientation_verified"))

print("\nPASS: Phase 18B schema audit completed")
