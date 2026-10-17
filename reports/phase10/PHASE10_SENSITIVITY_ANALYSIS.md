# Phase 10 — H3K27ac Enrichment Sensitivity Analysis

## Research objective

Evaluate the robustness of SREBF2 signature-associated H3K27ac
enrichment to gene biotype restriction, baseline expression
stratification, and flexible gene-length adjustment.

## Data

The analysis used 17,097 annotation-matched RNA-seq genes,
including 1,316 signature genes and 15,781 non-signature genes.

The H3K27ac data originated from an independent HepG2 CUT&Tag
experiment, while the transcriptional signature was derived
from HeLa SREBF2 perturbation RNA-seq.

## Methods

Three sensitivity-analysis families were evaluated:

1. Protein-coding-only logistic regression, adjusting for
   log10 gene length and log10(baseMean + 1).

2. Expression-quartile-specific logistic regression,
   adjusting for log10 gene length and gene biotype.

3. Logistic regression using natural cubic splines
   (3 degrees of freedom) for log10 gene length,
   additionally adjusting for baseline expression
   and gene biotype.

Each analysis evaluated promoter, gene-body, and any
H3K27ac peak overlap.

Across 18 exploratory tests, p-values were adjusted
using the Benjamini-Hochberg procedure.

## Results

### Protein-coding-only models

| Category | Adjusted OR | 95% CI | BH FDR |
|---|---:|---|---:|
| Promoter | 1.740 | 1.235–2.453 | 0.00216 |
| Gene body | 1.846 | 1.451–2.349 | 3.56e-06 |
| Any | 1.815 | 1.432–2.299 | 3.56e-06 |

### Expression-quartile models

| Quartile | Promoter OR | Gene-body OR | Any-overlap OR |
|---|---:|---:|---:|
| Q1 | 1.883 | 1.874 | 2.013 |
| Q2 | 2.019 | 1.650 | 1.647 |
| Q3 | 1.783 | 2.704 | 2.534 |
| Q4 | 0.801 | 2.299 | 2.319 |

All Q1–Q3 models were significant after BH correction.
None of the Q4 models passed the 0.05 FDR threshold.

### Flexible gene-length models

| Category | Adjusted OR | 95% CI | BH FDR |
|---|---:|---|---:|
| Promoter | 1.860 | 1.430–2.418 | 1.10e-05 |
| Gene body | 1.919 | 1.571–2.346 | 1.65e-09 |
| Any | 1.958 | 1.613–2.377 | 1.88e-10 |

All 18 models converged; 15 passed BH FDR < 0.05.

## Interpretation

The overall H3K27ac enrichment association is supported by
protein-coding-only and flexible gene-length analyses.

Expression-stratified results suggest that effect-size
estimates and statistical precision may vary across
baseline expression levels. Formal interaction testing
is required before concluding that enrichment differs
between expression strata.

## Limitations

- RNA-seq and CUT&Tag originated from different cell lines.
- The CUT&Tag assay does not measure direct SREBP2 binding.
- The analyses are observational and exploratory.
- Residual confounding and genomic dependencies remain possible.
- Expression-quartile analyses have different signature
  sample sizes and statistical precision.

## Reproducibility

Script:
scripts/cut_and_tag/phase10_sensitivity_analysis.R

Results:
results_public/cut_and_tag/phase10/SREBF2_H3K27ac_PHASE10_SENSITIVITY.tsv
