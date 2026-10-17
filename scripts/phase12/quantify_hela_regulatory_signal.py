#!/usr/bin/env python3
"""Quantify continuous HeLa ChIP-seq bigWig signals over GENCODE v47 genes."""

import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

BACKGROUND = Path(
    "results/cut_and_tag/integration/"
    "SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
)
ANNOTATION = Path(
    "data/processed/cut_and_tag/annotation/"
    "gencode_v47_gene_bodies.bed"
)
BIGWIG_DIR = Path("data/processed/phase12/bigwig")
OUTPUT = Path("results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv")
QC_OUTPUT = Path("reports/phase12/phase12_signal_quantification_qc.txt")

TRACKS = {
    "Input_CholDep": "GSM8258573_HInput.bw",
    "H3K27ac": "GSM8258574_HH3K27ac.bw",
    "SREBP2_FBS": "GSM8650746_FBS-IP_FKR0240286762-2A_hg38.bw",
    "SREBP2_CholDep": "GSM8650747_LPPS-IP_FKR0240286764-2A_hg38.bw",
}

PROMOTER_RADIUS = 2000


def measure(bw, chrom, start, end):
    """Return mean of defined values and defined-base fraction."""
    chrom_length = bw.chroms(chrom)

    if chrom_length is None:
        return np.nan, np.nan

    start = max(0, int(start))
    end = min(int(end), chrom_length)

    if end <= start:
        return np.nan, np.nan

    # exact=True avoids approximate bigWig zoom-level summaries.
    covered = bw.stats(
        chrom, start, end, type="coverage", exact=True
    )[0]

    if covered is None or not math.isfinite(covered):
        return np.nan, np.nan

    if covered <= 0:
        return np.nan, 0.0

    mean = bw.stats(
        chrom, start, end, type="mean", exact=True
    )[0]

    if mean is None or not math.isfinite(mean):
        return np.nan, float(covered)

    return float(mean), float(covered)


def main():
    background = pd.read_csv(
        BACKGROUND, sep="\t", dtype={"gene_id": str},
        low_memory=False
    )

    annotation = pd.read_csv(
        ANNOTATION,
        sep="\t",
        header=None,
        names=[
            "chrom", "start", "end", "gene_id",
            "annotation_symbol", "annotation_type", "strand"
        ],
        dtype={"gene_id": str},
    )

    if background["gene_id"].duplicated().any():
        raise ValueError("Duplicate gene IDs in RNA-seq background")

    if annotation["gene_id"].duplicated().any():
        raise ValueError("Duplicate gene IDs in GENCODE annotation")

    genes = background.merge(
        annotation,
        on="gene_id",
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    genes["gene_length"] = genes["end"] - genes["start"]

    tracks = {
        name: pyBigWig.open(str(BIGWIG_DIR / filename))
        for name, filename in TRACKS.items()
    }

    try:
        reference = tracks["H3K27ac"].chroms()

        for name, bw in tracks.items():
            if bw.chroms() != reference:
                raise ValueError(
                    f"Chromosome dictionary mismatch: {name}"
                )

        measurements = []

        for index, gene in enumerate(
            genes.itertuples(index=False), start=1
        ):
            result = {}

            if (
                pd.isna(gene.chrom)
                or pd.isna(gene.start)
                or pd.isna(gene.end)
                or gene.strand not in ("+", "-")
            ):
                for name in TRACKS:
                    result[f"mean_promoter_{name}"] = np.nan
                    result[f"promoter_defined_fraction_{name}"] = np.nan

                result["mean_genebody_H3K27ac"] = np.nan
                result["genebody_defined_fraction_H3K27ac"] = np.nan
                measurements.append(result)
                continue

            chrom = gene.chrom
            start = int(gene.start)
            end = int(gene.end)

            # BED coordinates: zero-based, half-open.
            # For minus-strand genes, the TSS base is end - 1.
            tss = start if gene.strand == "+" else end - 1

            promoter_start = tss - PROMOTER_RADIUS
            promoter_end = tss + PROMOTER_RADIUS

            for name, bw in tracks.items():
                mean, fraction = measure(
                    bw, chrom, promoter_start, promoter_end
                )
                result[f"mean_promoter_{name}"] = mean
                result[f"promoter_defined_fraction_{name}"] = fraction

            mean, fraction = measure(
                tracks["H3K27ac"], chrom, start, end
            )
            result["mean_genebody_H3K27ac"] = mean
            result["genebody_defined_fraction_H3K27ac"] = fraction

            measurements.append(result)

            if index % 2000 == 0:
                print(
                    f"Processed {index:,}/{len(genes):,} genes",
                    flush=True
                )

        measured = pd.DataFrame(measurements)
        output = pd.concat(
            [genes.reset_index(drop=True), measured],
            axis=1
        )

        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        output.to_csv(OUTPUT, sep="\t", index=False, na_rep="NA")

        signal_columns = [
            column for column in measured.columns
            if column.startswith("mean_")
        ]

        qc_lines = [
            "PHASE 12E — HELA REGULATORY SIGNAL QUANTIFICATION",
            f"Background rows: {len(background)}",
            f"Output rows: {len(output)}",
            f"Exact annotation matches: {(genes['_merge'] == 'both').sum()}",
            f"Unmatched annotation rows: {(genes['_merge'] != 'both').sum()}",
            f"Promoter radius: {PROMOTER_RADIUS} bp",
            "Promoter coordinates: BED zero-based, half-open",
            "Signal means: defined bases only; exact bigWig summaries",
            "",
            "SIGNAL COMPLETENESS",
        ]

        for column in signal_columns:
            qc_lines.append(
                f"{column}: "
                f"{output[column].notna().sum()} nonmissing; "
                f"{output[column].isna().sum()} missing"
            )

        QC_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        QC_OUTPUT.write_text("\n".join(qc_lines) + "\n")

        print("\n" + "\n".join(qc_lines))
        print(f"\nSaved: {OUTPUT}")
        print(f"Saved: {QC_OUTPUT}")

    finally:
        for bw in tracks.values():
            bw.close()


if __name__ == "__main__":
    main()
