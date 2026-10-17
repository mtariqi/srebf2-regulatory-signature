# SREBF2 Transcriptional Signature and H3K27ac Chromatin Context

**Analysis status:** Completed exploratory Phase 7–8 integration
**Date:** 17 October 2026
**Repository branch:** feature/phase7-cut-and-tag

## Research objective

Evaluate whether the frozen 1,326-gene SREBF2 perturbation-derived
transcriptional signature overlaps H3K27ac-marked genomic regions
more frequently than other RNA-seq-audited genes.

## Dataset provenance

- RNA-seq discovery: GSE267018, HeLa SREBF2 knockout.
- H3K27ac CUT&Tag: GSE228074, HepG2 control-shRNA condition.
- H3K27ac sample: SRR23949120.
- IgG control: SRR25870186.
- Reference annotation: GENCODE v47, GRCh38.

The CUT&Tag data do not represent SREBF2 perturbation or
SREBP2-specific chromatin occupancy.

## Methods

H3K27ac and IgG sequencing reads were subjected to quality control,
preprocessing, Bowtie2 alignment, duplicate processing, and MACS3
broad-peak calling.

Called peaks were associated with GENCODE v47 gene annotations.
Promoters were defined as transcription start sites +/- 2 kb.
Gene-body overlap was evaluated across annotated gene spans.

The frozen SREBF2 signature was integrated using versioned Ensembl
gene identifiers. Fisher's exact tests compared H3K27ac overlap
between selected and nonselected RNA-seq-audited genes.

P-values for the three overlap categories were corrected using
the Benjamini-Hochberg method.

## Integration quality control

| Metric | Count |
|---|---:|
| Frozen signature genes | 1,326 |
| Signature genes matched | 1,316 |
| Total RNA-seq-audited genes | 17,158 |
| Audited genes matched | 17,097 |
| Nonselected matched genes | 15,781 |
| Unmatched audited genes | 61 |

The ten unmatched signature genes included seven mitochondrial
genes and three lncRNAs.

## Background enrichment results

| Category | Selected overlap | Nonselected overlap | OR | 95% CI | BH FDR |
|---|---:|---:|---:|---|---:|
| Promoter | 75/1316 | 434/15781 | 2.137 | 1.638–2.757 | 4.268e-08 |
| Gene body | 239/1316 | 1591/15781 | 1.979 | 1.697–2.302 | 4.718e-17 |
| Any overlap | 249/1316 | 1655/15781 | 1.992 | 1.712–2.311 | 1.349e-17 |

## Direction-specific signature comparison

Genes increased after SREBF2 knockout exhibited higher odds of
H3K27ac overlap than genes decreased after knockout.

| Category | Increased-in-KO overlap | Decreased-in-KO overlap | OR | BH FDR |
|---|---:|---:|---:|---:|
| Promoter | 41/527 | 34/789 | 1.872 | 1.041e-02 |
| Gene body | 150/527 | 89/789 | 3.127 | 1.067e-14 |
| Any overlap | 155/527 | 94/789 | 3.078 | 1.067e-14 |

## Interpretation

The SREBF2 perturbation-derived transcriptional signature is
associated with increased prevalence of called H3K27ac peak
overlaps in an independent HepG2 chromatin dataset.

This is a cross-context association and does not demonstrate
direct SREBP2 binding, causal regulatory activity, or changes in
H3K27ac following SREBF2 perturbation.

## Limitations

- RNA-seq and CUT&Tag originate from different cell lines.
- The CUT&Tag dataset lacks a matched SREBF2 perturbation.
- The analyzed peak set lacks verified biological replication.
- The enrichment tests are unadjusted for gene length,
  baseline expression, and gene biotype.
- Larger genes may have more opportunities for peak overlap.
- Promoter and gene-body overlap categories are not independent.
- Absence of a called peak does not establish absence of signal.

## Next analysis

Phase 9 will assess whether the enrichment persists after
adjustment for gene length, RNA-seq expression, and gene biotype.

## Reproducibility

The statistical output tables are archived in:

`results_public/cut_and_tag/phase7_8/`

The analysis scripts are maintained under:

`scripts/cut_and_tag/`

Raw sequencing data and large intermediates are excluded from Git.
