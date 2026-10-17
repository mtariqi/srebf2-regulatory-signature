# Phase 9 — Covariate-Adjusted H3K27ac Enrichment

**Status:** Completed computational analysis
**Date:** 17 October 2026

## Objective

Determine whether the association between membership in the frozen
SREBF2 transcriptional signature and H3K27ac peak overlap persists
after adjustment for gene length, baseline RNA-seq expression,
and gene biotype.

## Methods

Binomial logistic regression was performed separately for
promoter, gene-body, and any H3K27ac peak overlap.

Predictors included:

- SREBF2 signature membership.
- Log10 gene length.
- Log10(baseMean + 1).
- Gene biotype (protein-coding, lncRNA, other).

The analysis included 17,097 annotation-matched genes:
1,316 signature and 15,781 non-signature genes.

Wald confidence intervals were calculated for the signature
coefficient. P-values across three overlap categories were
adjusted using the Benjamini-Hochberg procedure.

## Results

| Category | Adjusted OR | 95% CI | BH FDR |
|---|---:|---|---:|
| Promoter | 1.938 | 1.492–2.516 | 6.898e-07 |
| Gene body | 1.926 | 1.580–2.347 | 1.293e-10 |
| Any overlap | 1.971 | 1.629–2.385 | 9.329e-12 |

All three logistic regression models converged.

## Interpretation

Signature genes exhibited higher odds of H3K27ac peak overlap
than non-signature genes after adjustment for the specified
covariates.

The results support an association between the SREBF2
perturbation-derived transcriptional signature and H3K27ac-marked
chromatin in an independent dataset.

## Limitations

RNA-seq originated from HeLa cells, whereas H3K27ac CUT&Tag
originated from HepG2 cells.

The CUT&Tag data do not represent SREBF2 perturbation or
direct SREBP2 occupancy.

Adjustment for gene length, expression, and biotype does not
eliminate all potential confounding. The models do not establish
causal regulatory effects.

## Reproducibility

Script:
scripts/cut_and_tag/adjusted_h3k27ac_enrichment.R

Statistical results:
results_public/cut_and_tag/phase9/SREBF2_H3K27ac_ADJUSTED_ENRICHMENT.tsv

Model QC:
reports/phase9/adjusted_model_QC.txt
