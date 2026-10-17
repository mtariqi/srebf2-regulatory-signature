# Phase 11 — Expression Interaction and Direction-Specific Enrichment

## Research objective

Determine whether H3K27ac enrichment associated with the
SREBF2 transcriptional signature varies by baseline expression
or by transcriptional response direction.

## Dataset

17,097 annotation-matched genes:

- 15,781 non-signature genes
- 789 decreased-in-KO signature genes
- 527 increased-in-KO signature genes

## Methods

### Expression interaction

Nested binomial logistic regression models with and without
signature-by-expression-quartile interaction were compared
using likelihood-ratio tests.

Both models adjusted for gene length and gene biotype.

### Direction-specific enrichment

Multivariable logistic regression compared each signature
direction with non-signature genes.

Covariates:

- log10 gene length
- log10(baseMean + 1)
- Gene biotype

### Direct direction contrasts

Increased-in-KO and decreased-in-KO coefficients were
compared using the fitted model covariance matrix.

Wald confidence intervals and two-sided p-values were
calculated for the adjusted odds-ratio contrasts.

Benjamini-Hochberg adjustment was applied separately to
the three interaction tests, six direction-specific tests,
and three direct contrasts.

## Results

### Expression interaction

| Category | LR chi-square | df | BH FDR |
|---|---:|---:|---:|
| Promoter | 0.897 | 3 | 0.826 |
| Gene body | 3.983 | 3 | 0.487 |
| Any overlap | 3.469 | 3 | 0.487 |

No statistically significant interaction was detected.

### Direction-specific enrichment

| Category | Direction | Adjusted OR | BH FDR |
|---|---|---:|---:|
| Promoter | Decreased in KO | 1.946 | 3.70e-04 |
| Promoter | Increased in KO | 1.931 | 2.29e-04 |
| Gene body | Decreased in KO | 2.354 | 3.67e-08 |
| Gene body | Increased in KO | 1.684 | 6.86e-05 |
| Any overlap | Decreased in KO | 2.324 | 2.46e-08 |
| Any overlap | Increased in KO | 1.755 | 1.15e-05 |

All six direction-specific associations were significant
after BH correction.

### Direct increased-versus-decreased contrasts

| Category | Adjusted OR | 95% CI | BH FDR |
|---|---:|---|---:|
| Promoter | 0.992 | 0.612–1.609 | 0.975 |
| Gene body | 0.715 | 0.492–1.040 | 0.190 |
| Any overlap | 0.755 | 0.527–1.083 | 0.190 |

None of the direct direction contrasts passed BH FDR 0.05.

## Scientific interpretation

The SREBF2 transcriptional signature exhibits H3K27ac
enrichment relative to non-signature genes for both
transcriptional response directions.

There was no statistically significant evidence of
expression-quartile interaction or a difference between
the two signature directions.

Nonsignificant tests do not establish equivalence or
prove the absence of biological differences.

## Limitations

The RNA-seq signature originated from HeLa SREBF2
perturbation data, whereas the H3K27ac CUT&Tag dataset
originated from HepG2 cells.

The CUT&Tag experiment does not directly measure SREBP2
binding or SREBF2 perturbation effects.

The analyses remain exploratory and observational.
Residual confounding and genomic dependencies may remain.

## Reproducibility

Scripts:

- scripts/cut_and_tag/phase11_interaction_direction.R
- scripts/cut_and_tag/phase11_direction_contrasts.R

Results:

- results_public/cut_and_tag/phase11/SREBF2_H3K27ac_EXPRESSION_INTERACTION.tsv
- results_public/cut_and_tag/phase11/SREBF2_H3K27ac_DIRECTION_ADJUSTED.tsv
- results_public/cut_and_tag/phase11/SREBF2_H3K27ac_DIRECTION_CONTRASTS.tsv
