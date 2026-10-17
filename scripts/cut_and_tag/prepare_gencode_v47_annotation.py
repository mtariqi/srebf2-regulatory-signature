#!/usr/bin/env python3

import re
from pathlib import Path
from collections import Counter

gtf = Path(
    "data/raw/reference/gencode_v47/"
    "gencode.v47.primary_assembly.annotation.gtf"
)

outdir = Path("data/processed/cut_and_tag/annotation")
outdir.mkdir(parents=True, exist_ok=True)

promoter_file = outdir / "gencode_v47_promoters_2kb.bed"
gene_file = outdir / "gencode_v47_gene_bodies.bed"

standard = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX", "chrY"
}

def get_attribute(attributes, key):
    match = re.search(
        rf'(?:^|;\s*){re.escape(key)} "([^"]+)"',
        attributes
    )
    return match.group(1) if match else "NA"

counts = Counter()

with (
    gtf.open() as source,
    promoter_file.open("w") as promoters,
    gene_file.open("w") as genes
):
    for line in source:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9 or fields[2] != "gene":
            continue

        chrom = fields[0]

        if chrom not in standard:
            counts["excluded_nonstandard_genes"] += 1
            continue

        start_1 = int(fields[3])
        end_1 = int(fields[4])
        strand = fields[6]

        if strand not in ("+", "-"):
            counts["excluded_unknown_strand"] += 1
            continue

        attributes = fields[8]

        gene_id = get_attribute(attributes, "gene_id")
        gene_name = get_attribute(attributes, "gene_name")
        gene_type = get_attribute(attributes, "gene_type")

        # GTF: 1-based inclusive.
        # BED: 0-based half-open.
        gene_start = start_1 - 1
        gene_end = end_1

        # 0-based coordinate of the transcription start site.
        tss = start_1 - 1 if strand == "+" else end_1 - 1

        promoter_start = max(0, tss - 2000)
        promoter_end = tss + 2001

        metadata = [
            gene_id,
            gene_name,
            gene_type,
            strand
        ]

        promoters.write(
            "\t".join(
                map(str, [
                    chrom,
                    promoter_start,
                    promoter_end,
                    *metadata
                ])
            ) + "\n"
        )

        genes.write(
            "\t".join(
                map(str, [
                    chrom,
                    gene_start,
                    gene_end,
                    *metadata
                ])
            ) + "\n"
        )

        counts["retained_genes"] += 1

print("GENCODE v47 annotation preparation completed")
print(f"Retained genes: {counts['retained_genes']:,}")
print(
    "Excluded nonstandard-chromosome genes:",
    f"{counts['excluded_nonstandard_genes']:,}"
)
print(
    "Excluded unknown-strand genes:",
    f"{counts['excluded_unknown_strand']:,}"
)
print(f"Promoter BED: {promoter_file}")
print(f"Gene-body BED: {gene_file}")
