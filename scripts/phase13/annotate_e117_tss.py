#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase13/annotate_e117_tss.py
# Purpose     : Annotate frozen SREBF2 gene universe with Roadmap E117 TSS states
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Annotate frozen SREBF2 gene universe with Roadmap E117 TSS states."""

import csv
import gzip
import hashlib
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path

MATRIX = Path(
    "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
)
SEGMENTATION = Path(
    "metadata/phase13/reference/"
    "E117_15_coreMarks_hg38lift_mnemonics.bed.gz"
)
OUTPUT = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)

EXPECTED_SHA256 = (
    "f1b30d0fa472c61642b35e2819d130a"
    "ec0b0359d9bf6abdbdf66418f18271508"
)

EXPECTED_GROUPS = {
    "decreased_in_KO": 789,
    "increased_in_KO": 527,
    "background": 15781,
}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_coordinate(value):
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"Noninteger coordinate: {value}")
    return int(number)

def main():
    if OUTPUT.exists():
        raise FileExistsError(
            f"Output already exists; refusing overwrite: {OUTPUT}"
        )

    actual_hash = sha256_file(SEGMENTATION)
    if actual_hash != EXPECTED_SHA256:
        raise ValueError("Reference SHA-256 mismatch")

    intervals = defaultdict(list)

    with gzip.open(SEGMENTATION, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, state = line.rstrip().split("\t")[:4]
            intervals[chrom].append((int(start), int(end), state))

    indexes = {}

    for chrom, records in intervals.items():
        records.sort(key=lambda x: x[0])

        for i, (start, end, _) in enumerate(records):
            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval on {chrom}")
            if i and start < records[i - 1][1]:
                raise ValueError(f"Overlapping intervals on {chrom}")

        starts = [record[0] for record in records]
        indexes[chrom] = (starts, records)

    def assign_state(chrom, tss):
        if chrom not in indexes:
            return "UNANNOTATED"

        starts, records = indexes[chrom]
        i = bisect_right(starts, tss) - 1

        if i >= 0:
            start, end, state = records[i]
            if start <= tss < end:
                return state

        return "UNANNOTATED"

    rows = []
    groups = Counter()
    states = Counter()

    with open(MATRIX, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for gene in reader:
            if gene["strand"] not in ("+", "-"):
                continue

            start = parse_coordinate(gene["start"])
            end = parse_coordinate(gene["end"])

            if start < 0 or end <= start:
                raise ValueError(
                    f"Invalid gene coordinates: {gene['gene_id']}"
                )

            tss = start if gene["strand"] == "+" else end - 1

            selected = gene["signature_selected"] == "True"

            if selected:
                group = gene["direction"]
                if group not in (
                    "decreased_in_KO",
                    "increased_in_KO",
                ):
                    raise ValueError(
                        f"Unexpected signature direction: {group}"
                    )
            else:
                group = "background"

            state = assign_state(gene["chrom"], tss)

            rows.append({
                "gene_id": gene["gene_id"],
                "gene_symbol": gene["gene_symbol"],
                "gene_type": gene["gene_type"],
                "chrom": gene["chrom"],
                "strand": gene["strand"],
                "tss_0based": tss,
                "chromhmm_state": state,
                "signature_group": group,
                "baseMean": gene["baseMean"],
                "gene_length": gene["gene_length"],
            })

            groups[group] += 1
            states[state] += 1

    if dict(groups) != EXPECTED_GROUPS:
        raise ValueError(
            f"Unexpected group counts: {dict(groups)}"
        )

    fieldnames = [
        "gene_id", "gene_symbol", "gene_type",
        "chrom", "strand", "tss_0based",
        "chromhmm_state", "signature_group",
        "baseMean", "gene_length",
    ]

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with open(OUTPUT, "x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print("=== Phase 13C: E117 TSS annotation ===")
    print("Status: PASS")
    print("Total annotated gene records:", len(rows))
    print("Reference SHA-256:", actual_hash)

    print("\nSignature groups:")
    for group, count in sorted(groups.items()):
        print(group, count)

    print("\nChromatin-state assignments:")
    for state, count in sorted(states.items()):
        print(state, count)

    print("\nOutput:", OUTPUT)

if __name__ == "__main__":
    main()
