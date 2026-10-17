# Bowtie2 Dovetail Parameter Validation

## Dataset
- SRA accession: SRR23949120
- Assay: H3K27ac CUT&Tag
- Experimental condition: HepG2 control shRNA
- Reference: GRCh38
- Pilot input: First 100,000 paired-end read pairs

## Experimental comparison
Identical input reads were aligned using Bowtie2
with and without the --dovetail option.

Common parameters:
--very-sensitive --no-mixed --no-discordant -X 2000 -p 2

## Results

| Metric | Dovetail | No dovetail |
|---|---:|---:|
| Complete primary pairs | 94094 | 92396 |
| Both mates MAPQ >= 30 | 84532 | 83950 |
| Incomplete/invalid groups | 0 | 0 |
| High-MAPQ pair fraction | 89.84% | 90.86% |
| Positive-TLEN records below 200 bp | 9636 | 8782 |
| Minimum observed positive TLEN | 134 bp | 150 bp |

## Interpretation

Allowing dovetailing increased recovery of complete primary
pairs and short-fragment records. The absolute number of
pairs with MAPQ >= 30 increased, although their proportion
decreased slightly.

Dovetailing is provisionally selected for the final
alignment workflow.

## Limitations

This pilot represents only the first 100,000 read pairs
from one sequencing library.

Further assessment should include duplicate rates,
fragment-size distributions, and enrichment relative
to matched IgG controls.

These H3K27ac CUT&Tag data investigate KDM1A-associated
chromatin regulation, not direct SREBF2 occupancy.
