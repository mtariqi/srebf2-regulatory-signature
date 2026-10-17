#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

signature_file = Path(
    "results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv"
)

h3k27ac_file = Path(
    "results/cut_and_tag/annotation/H3K27ac_gene_level_summary.tsv"
)

output = Path(
    "results/cut_and_tag/integration/"
    "SREBF2_SIGNATURE_H3K27ac_INTEGRATED.tsv"
)

report = Path(
    "reports/phase7/SREBF2_H3K27ac_integration_summary.txt"
)

output.parent.mkdir(parents=True, exist_ok=True)
report.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stable_id(gene_id):
    return gene_id.split(".")[0]


signature = read_tsv(signature_file)
chromatin = read_tsv(h3k27ac_file)

exact = {}
stable = defaultdict(list)

for row in chromatin:
    gene_id = row["gene_id"]

    if gene_id in exact:
        raise ValueError(f"Duplicate annotation gene ID: {gene_id}")

    exact[gene_id] = row
    stable[stable_id(gene_id)].append(row)

if len({r["gene_id"] for r in signature}) != len(signature):
    raise ValueError("Duplicate gene IDs in signature")

columns = list(signature[0].keys()) + [
    "h3k27ac_gene_id",
    "h3k27ac_gene_name",
    "h3k27ac_chrom",
    "h3k27ac_strand",
    "promoter_peak_count",
    "gene_body_peak_count",
    "any_peak_count",
    "h3k27ac_match_status",
]

counts = Counter()
direction_counts = defaultdict(Counter)
integrated = []

for row in signature:
    gene_id = row["gene_id"]

    if gene_id in exact:
        match = exact[gene_id]
        status = "exact"

    else:
        candidates = stable.get(stable_id(gene_id), [])

        if len(candidates) == 1:
            match = candidates[0]
            status = "stable_id"

        elif len(candidates) > 1:
            match = None
            status = "ambiguous"

        else:
            match = None
            status = "unmatched"

    result = dict(row)

    if match is not None:
        result.update({
            "h3k27ac_gene_id": match["gene_id"],
            "h3k27ac_gene_name": match["gene_name"],
            "h3k27ac_chrom": match["chrom"],
            "h3k27ac_strand": match["strand"],
            "promoter_peak_count": match["promoter_peak_count"],
            "gene_body_peak_count": match["gene_body_peak_count"],
            "any_peak_count": match["any_peak_count"],
        })
    else:
        result.update({
            "h3k27ac_gene_id": "",
            "h3k27ac_gene_name": "",
            "h3k27ac_chrom": "",
            "h3k27ac_strand": "",
            "promoter_peak_count": "",
            "gene_body_peak_count": "",
            "any_peak_count": "",
        })

    result["h3k27ac_match_status"] = status
    integrated.append(result)

    counts[status] += 1

    if match is not None:
        direction = row["direction"]

        direction_counts[direction]["total_matched"] += 1

        if int(match["promoter_peak_count"]) > 0:
            direction_counts[direction]["promoter"] += 1

        if int(match["gene_body_peak_count"]) > 0:
            direction_counts[direction]["gene_body"] += 1

        if int(match["any_peak_count"]) > 0:
            direction_counts[direction]["any"] += 1


with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(integrated)


lines = [
    "SREBF2 SIGNATURE × H3K27ac INTEGRATION",
    "=====================================",
    f"Signature genes: {len(signature):,}",
    f"Exact Ensembl matches: {counts['exact']:,}",
    f"Stable-ID matches: {counts['stable_id']:,}",
    f"Ambiguous matches: {counts['ambiguous']:,}",
    f"Unmatched genes: {counts['unmatched']:,}",
    "",
    "H3K27ac OVERLAPS BY SIGNATURE DIRECTION",
]

for direction in sorted(direction_counts):
    c = direction_counts[direction]

    lines.extend([
        "",
        f"Direction: {direction}",
        f"Matched genes: {c['total_matched']:,}",
        f"Genes with promoter peaks: {c['promoter']:,}",
        f"Genes with gene-body peaks: {c['gene_body']:,}",
        f"Genes with any peaks: {c['any']:,}",
    ])

lines.extend([
    "",
    f"Integrated table: {output}",
])

report.write_text("\n".join(lines) + "\n")

print("\n".join(lines))

if len(integrated) != 1326:
    raise SystemExit(
        f"ERROR: Expected 1,326 signature genes, "
        f"found {len(integrated)}"
    )
