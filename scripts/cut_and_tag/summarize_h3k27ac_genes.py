#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("results/cut_and_tag/annotation")
annot = Path("data/processed/cut_and_tag/annotation")

promoter_overlaps = base / "H3K27ac_promoter_overlaps.tsv"
body_overlaps = base / "H3K27ac_gene_body_overlaps.tsv"
gene_bed = annot / "gencode_v47_gene_bodies.bed"

output = base / "H3K27ac_gene_level_summary.tsv"

genes = {}
counts = defaultdict(
    lambda: {
        "promoter": set(),
        "body": set()
    }
)

# Load the complete gene reference.
with gene_bed.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")

        chrom, start, end, gene_id, name, biotype, strand = fields

        genes[gene_id] = {
            "gene_name": name,
            "gene_type": biotype,
            "chrom": chrom,
            "strand": strand
        }

def load_overlaps(path, category):
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            # MACS3 broadPeak contains 9 columns.
            # BEDTools -wa -wb appends 7 annotation columns.
            if len(fields) != 16:
                raise ValueError(
                    f"Expected 16 columns in {path}; found {len(fields)}"
                )

            peak_chrom = fields[0]
            peak_start = fields[1]
            peak_end = fields[2]

            gene_id = fields[12]

            if gene_id not in genes:
                raise ValueError(f"Unknown gene ID: {gene_id}")

            peak_key = (peak_chrom, peak_start, peak_end)

            counts[gene_id][category].add(peak_key)

load_overlaps(promoter_overlaps, "promoter")
load_overlaps(body_overlaps, "body")

columns = [
    "gene_id",
    "gene_name",
    "gene_type",
    "chrom",
    "strand",
    "promoter_peak_count",
    "gene_body_peak_count",
    "any_peak_count"
]

with output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(columns)

    for gene_id in sorted(genes):
        metadata = genes[gene_id]

        promoter = counts[gene_id]["promoter"]
        body = counts[gene_id]["body"]

        writer.writerow([
            gene_id,
            metadata["gene_name"],
            metadata["gene_type"],
            metadata["chrom"],
            metadata["strand"],
            len(promoter),
            len(body),
            len(promoter | body)
        ])

associated = sum(
    bool(v["promoter"] or v["body"])
    for v in counts.values()
)

promoter_genes = sum(
    bool(v["promoter"])
    for v in counts.values()
)

body_genes = sum(
    bool(v["body"])
    for v in counts.values()
)

print("H3K27ac gene-level summary completed")
print(f"Total annotated genes: {len(genes):,}")
print(f"Genes with any overlap: {associated:,}")
print(f"Genes with promoter overlap: {promoter_genes:,}")
print(f"Genes with gene-body overlap: {body_genes:,}")
print(f"Output: {output}")
