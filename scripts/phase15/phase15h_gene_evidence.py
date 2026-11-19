#!/usr/bin/env python3
"""Build a gene-level evidence table from frozen SREBF2 analyses."""

from pathlib import Path
import pandas as pd

BASE = Path("results/phase15")
ENRICH = BASE / "phase15e_enrichment"
OUT = BASE / "phase15h_gene_evidence.tsv"

if OUT.exists():
    raise SystemExit(f"Refusing to overwrite: {OUT}")

pathways = {
    "GO_sterol": ("GO_BP_GSEA", "GO:0016126"),
    "GO_cholesterol": ("GO_BP_GSEA", "GO:0006695"),
    "Reactome_cholesterol": (
        "Reactome_GSEA", "R-HSA-191273"
    ),
    "Reactome_SREBF": (
        "Reactome_GSEA", "R-HSA-2426168"
    ),
    "KEGG_steroid": ("KEGG_GSEA", "hsa00100"),
    "KEGG_terpenoid": ("KEGG_GSEA", "hsa00900"),
    "KEGG_oxidative": ("KEGG_GSEA", "hsa00190"),
}

mapping = pd.read_csv(
    ENRICH / "ensembl_entrez_mapping.tsv",
    sep="\t",
    dtype=str,
)

de = pd.read_csv(
    "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv",
    sep="\t",
)

de["ENSEMBL"] = (
    de["gene_id"].astype(str)
    .str.replace(r"\.\d+$", "", regex=True)
)

if de["ENSEMBL"].duplicated().any():
    raise ValueError("Duplicate Ensembl IDs in Phase 12")

if mapping["ENTREZID"].duplicated().any():
    raise ValueError("Duplicate Entrez IDs in mapping")

cols = [
    "ENSEMBL",
    "gene_symbol",
    "log2FC_primary",
    "padj_primary",
    "log2FC_fbs",
    "padj_fbs",
    "log2FC_sensitivity",
    "padj_sensitivity",
]

genes = mapping.merge(
    de[cols],
    on="ENSEMBL",
    how="left",
    validate="one_to_one",
)

if genes["gene_symbol"].isna().any():
    raise ValueError("Missing Phase 12 gene mappings")

membership = {}

for label, (analysis, pathway_id) in pathways.items():
    table = pd.read_csv(
        ENRICH / f"{analysis}.tsv",
        sep="\t",
        dtype={"ID": str},
    )

    row = table.loc[table["ID"] == pathway_id]

    if len(row) != 1:
        raise ValueError(f"Pathway missing: {label}")

    membership[label] = set(
        str(row.iloc[0]["core_enrichment"]).split("/")
    )

    genes[label] = genes["ENTREZID"].isin(
        membership[label]
    )

genes["leading_edge_pathway_count"] = (
    genes[list(pathways)].sum(axis=1)
)

genes = genes.loc[
    genes["leading_edge_pathway_count"] > 0
].copy()

genes = genes.sort_values(
    ["leading_edge_pathway_count", "log2FC_primary"],
    ascending=[False, True],
)

genes.to_csv(
    OUT,
    sep="\t",
    index=False,
    mode="x",
)

print("PASS: Phase 15H gene evidence generated")
print("Output:", OUT)
print("Unique genes:", len(genes))

core = (
    membership["GO_sterol"]
    & membership["Reactome_cholesterol"]
    & membership["KEGG_steroid"]
)

print("\nShared three-pathway core:", len(core))

print(
    genes.loc[
        genes["ENTREZID"].isin(core),
        [
            "gene_symbol",
            "log2FC_primary",
            "padj_primary",
            "leading_edge_pathway_count",
        ],
    ]
    .sort_values("log2FC_primary")
    .to_string(index=False)
)

print("\nAll Phase 12 and Phase 15E inputs unchanged.")
