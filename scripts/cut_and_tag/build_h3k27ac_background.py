#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

audit_file = Path(
    "results/signature/discovery_v1/"
    "SIGNATURE_SELECTION_AUDIT.tsv"
)

annotation_file = Path(
    "results/cut_and_tag/annotation/"
    "H3K27ac_gene_level_summary.tsv"
)

output_file = Path(
    "results/cut_and_tag/integration/"
    "SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
)

report_file = Path(
    "reports/phase7/"
    "SREBF2_H3K27ac_background_QC.txt"
)

output_file.parent.mkdir(parents=True, exist_ok=True)
report_file.parent.mkdir(parents=True, exist_ok=True)


def load_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stable_id(value):
    return value.split(".")[0]


audit = load_tsv(audit_file)
annotation = load_tsv(annotation_file)

exact = {}
stable = defaultdict(list)

for gene in annotation:
    gene_id = gene["gene_id"]

    if gene_id in exact:
        raise ValueError(f"Duplicate annotation ID: {gene_id}")

    exact[gene_id] = gene
    stable[stable_id(gene_id)].append(gene)

if len({r["gene_id"] for r in audit}) != len(audit):
    raise ValueError("Duplicate gene IDs in RNA-seq audit")

extra_columns = [
    "h3k27ac_gene_id",
    "h3k27ac_chrom",
    "h3k27ac_strand",
    "promoter_peak_count",
    "gene_body_peak_count",
    "any_peak_count",
    "h3k27ac_match_status",
]

fields = list(audit[0]) + extra_columns

qc = Counter()
output_rows = []

for gene in audit:
    gene_id = gene["gene_id"]

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

    row = dict(gene)

    if match:
        row.update({
            "h3k27ac_gene_id": match["gene_id"],
            "h3k27ac_chrom": match["chrom"],
            "h3k27ac_strand": match["strand"],
            "promoter_peak_count": match["promoter_peak_count"],
            "gene_body_peak_count": match["gene_body_peak_count"],
            "any_peak_count": match["any_peak_count"],
        })
    else:
        row.update({
            key: "" for key in extra_columns[:-1]
        })

    row["h3k27ac_match_status"] = status

    qc[status] += 1
    qc["total"] += 1

    selected = gene["signature_selected"] == "True"
    qc["selected" if selected else "nonselected"] += 1

    if match:
        group = "selected" if selected else "nonselected"
        qc[f"{group}_matched"] += 1

        for category in (
            "promoter_peak_count",
            "gene_body_peak_count",
            "any_peak_count"
        ):
            if int(match[category]) > 0:
                qc[f"{group}_{category}_positive"] += 1

    output_rows.append(row)


with output_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output_rows)

lines = [
    "SREBF2 RNA-SEQ × H3K27ac BACKGROUND QC",
    "====================================",
    f"Total audited genes: {qc['total']:,}",
    f"Selected signature genes: {qc['selected']:,}",
    f"Nonselected genes: {qc['nonselected']:,}",
    "",
    f"Exact matches: {qc['exact']:,}",
    f"Stable-ID matches: {qc['stable_id']:,}",
    f"Ambiguous matches: {qc['ambiguous']:,}",
    f"Unmatched genes: {qc['unmatched']:,}",
]

for group in ("selected", "nonselected"):
    lines.extend([
        "",
        f"{group.upper()}",
        f"Matched genes: {qc[f'{group}_matched']:,}",
        "Promoter-positive: "
        f"{qc[f'{group}_promoter_peak_count_positive']:,}",
        "Gene-body-positive: "
        f"{qc[f'{group}_gene_body_peak_count_positive']:,}",
        "Any-overlap-positive: "
        f"{qc[f'{group}_any_peak_count_positive']:,}",
    ])

report_file.write_text("\n".join(lines) + "\n")

print("\n".join(lines))

assert qc["total"] == 17158
assert qc["selected"] == 1326
assert qc["nonselected"] == 15832
