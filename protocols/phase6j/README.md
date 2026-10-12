# Phase 6J — SREBF2 Experimental Validation Protocol

**Status:** Prospective experimental plan; no wet-lab validation has yet been performed.

## 1. Scientific rationale

Independent HAP1 RNA-seq analysis (GSE147770) identified genotype-by-serum interactions among a frozen SREBF2-associated signature. Of 181 genes showing discovery-consistent +FBS effects, opposite −FBS effects, and FDR-significant interactions, 13 overlapped the Hallmark oxidative phosphorylation gene set. Phase 6I prioritized five candidates using RNA interaction strength and reversal magnitude. Promoter-associated SREBP2 ChIP-seq coverage (GSE282800) is exploratory and did not demonstrate consistently stronger signals after matched-background comparison. Cholesterol depletion in the ChIP-seq dataset is **not** the same treatment as serum starvation in HAP1.

## 2. Primary aims

1. Replicate genotype-by-serum interaction effects for ETFB, TIMM10, MRPS12, NDUFB7, and ATP5MC3 by RT-qPCR in independently prepared HAP1 cultures.
2. Determine whether SREBF2 specifically occupies candidate regulatory regions using an independently controlled occupancy assay.
3. Determine whether genotype-by-serum interaction is accompanied by a change in mitochondrial respiration, without assuming that transcript changes predict function.
4. Evaluate genotype specificity by complementation/rescue and, if feasible, independent SREBF2 perturbation.

## 3. Prespecified hypotheses

- **H1:** Each candidate has a nonzero genotype-by-serum interaction on normalized log-expression in HAP1 cells.
- **H2:** The KO-vs-WT effect is lower under +FBS than under −FBS, consistent with the RNA-seq interaction direction.
- **H3 (exploratory):** Some candidate promoters exhibit specific SREBF2 occupancy relative to validated positive and negative controls.
- **H4 (exploratory):** Respiration-related phenotypes show genotype-by-serum interactions and may be partially rescued by SREBF2 restoration.

## 4. Primary RT-qPCR design

| Group | Genotype | Serum | Exposure | Independent biological replicates |
|---|---|---|---|---|
| A | WT HAP1 | +FBS | matched harvest | 3 pilot; target >=5–6 for confirmation after power assessment |
| B | SREBF2 KO HAP1 | +FBS | matched harvest | same |
| C | WT HAP1 | −FBS | 4 hours, confirm against original GSE147770 methods | same |
| D | SREBF2 KO HAP1 | −FBS | 4 hours, confirm against original GSE147770 methods | same |

Biological replicates must be independent culture preparations, ideally across days/batches. Balance batch across the four groups, randomize processing order, and blind qPCR plate annotation where practical. Technical wells are not independent biological replicates. Include no-template and no-RT controls. Confirm mycoplasma-negative status, cell identity, passage range, and comparable confluence/viability. Use an appropriate matched control medium for the serum-withdrawal comparison; document medium, serum lot, timing, and density. Verify knockout genotype and SREBF2 protein status; residual RNA alone does not demonstrate functional protein.

## 5. Prespecified RNA candidates

| Gene | KO-vs-WT log2FC +FBS | KO-vs-WT log2FC −FBS | Δ(−FBS − +FBS) | DESeq2 interaction FDR |
|---|---:|---:|---:|---:|
| ETFB | −0.197352 | +0.594576 | +0.791928 | 2.08188e−10 |
| TIMM10 | −0.146681 | +0.573045 | +0.719725 | 3.10915e−08 |
| MRPS12 | −0.156933 | +0.445955 | +0.602888 | 3.91058e−08 |
| NDUFB7 | −0.088947 | +0.526975 | +0.615921 | 5.59509e−08 |
| ATP5MC3 | −0.195720 | +0.223852 | +0.419572 | 2.27948e−08 |

The values are computational estimates, not observed experimental validation results. The small +FBS effects make detection of condition-specific changes challenging even when interactions are strong.

## 6. RT-qPCR assay qualification and processing

1. Select exon-spanning assays where feasible; confirm transcript isoform coverage and specificity against the actual HAP1 genotype. Record primer/probe identifiers and sequences only after experimental verification.
2. Evaluate amplification specificity (melt curve for dye-based assays), efficiency across a dilution series, dynamic range, and technical variability. Follow MIQE principles.
3. Test multiple candidate reference genes across **all four groups**; select stable references empirically. Do not assume GAPDH or ACTB is stable during serum deprivation.
4. Use consistent RNA isolation, DNase treatment, RNA integrity/quality checks, input quantity, and reverse-transcription workflow. Randomize plate positions and include interplate controls if needed.
5. Predefine exclusion criteria for failed wells and samples without using the observed interaction effect to determine exclusions.
6. Preserve raw Ct values, amplification QC, sample metadata, assay efficiencies, and analysis scripts.

## 7. Primary statistical analysis

Fit each gene on normalized log2-expression, or an appropriately specified ΔCt scale, using:

`expression ~ genotype * serum + batch` (omit batch only if no independent batch effect is estimable).

Reference levels: WT and +FBS. The interaction coefficient is the KO effect in −FBS minus the KO effect in +FBS. Report its estimate, 95% confidence interval, raw p-value, and BH-adjusted p-value across the five prespecified candidate genes. For ΔCt, remember that increasing ΔCt indicates lower expression; reverse the sign when interpreting log2 expression. Treat qPCR technical replicates as nested measurements or summarize per biological replicate after QC. Use robust diagnostics and account for culture batch when justified. Do not claim confirmation based solely on matching fold-change signs; assess effect uncertainty. Plan confirmatory sample size using pilot residual variance and a prespecified minimum meaningful interaction effect, not only the very small RNA-seq p-values.

## 8. Binding occupancy validation

Prioritize ETFB and NDUFB7 for a small initial feasibility panel, with TIMM10, MRPS12, and ATP5MC3 following if controls and assay performance are adequate. Obtain the correct HAP1 genomic coordinates and promoter/transcription-start-site annotation before designing assays. Use SREBF2 CUT&RUN or ChIP-qPCR with antibody validation, input/background controls, known SREBF2-positive loci (e.g., a validated sterol-responsive locus in the same cell context), and negative genomic loci. Include WT and KO as specificity controls and compare matched +FBS/−FBS conditions. Avoid interpreting the separate HeLa cholesterol-depletion bigWig experiment as direct HAP1 serum-condition occupancy evidence. Predefine enrichment normalization, technical QC, replicate structure, and multiple-testing treatment.

## 9. Functional mitochondrial validation

Measure basal oxygen consumption, ATP-linked respiration, maximal respiration and spare respiratory capacity using a validated mitochondrial stress-test platform, with assay settings optimized for HAP1. Include cell number/protein normalization and viability, growth, and mitochondrial-content controls. Compare the same four genotype-by-serum groups and test the interaction; consider ATP measurements and mitochondrial membrane potential as secondary endpoints. Changes in respiratory transcripts alone cannot establish functional changes.

## 10. Rescue and specificity

Where feasible, complement the SREBF2 KO with physiologically controlled SREBF2 expression, with matched empty-vector controls and documented expression levels. Distinguish full-length SREBF2 from constitutively active nuclear forms; these constructs test different biological mechanisms. Consider an independent KO clone or orthogonal knockdown to assess clonal/off-target effects. Rescue must be evaluated for interaction effects, not merely one condition's expression.

## 11. Decision gates

- **Gate A — Technical:** Confirm cell-line identity, SREBF2 perturbation, stable references, and qPCR assay performance.
- **Gate B — RNA:** Quantify interaction estimates for five candidates; BH correction and confidence intervals; compare directions with Phase 6D.
- **Gate C — Occupancy:** Require specific, control-supported binding signal before using the phrase “direct SREBF2 target.”
- **Gate D — Function:** Require appropriately normalized respiration measurements and replicate-supported interaction effects before claiming altered OXPHOS function.
- **Gate E — Specificity:** Assess rescue and independent perturbation where possible.

## 12. Reproducibility package

Suggested repository destination: `protocols/phase6j/`.

Track `experimental_design.tsv`, `sample_manifest.tsv`, `qpcr_assay_validation.tsv`, `raw_ct.tsv`, `qpcr_qc.tsv`, `normalized_expression.tsv`, `statistical_analysis.R`, `binding_assay_manifest.tsv`, `mitochondrial_assay_manifest.tsv`, `README.md`, and an analysis environment/session record. Populate experimental tables only after observations exist. Do not fabricate primer sequences, measured Ct values, binding occupancy, or mitochondrial assay outputs.

## 13. Interpretation boundary

The computational evidence establishes serum-dependent transcriptional interaction and exploratory mitochondrial pathway enrichment. It does **not** establish direct SREBF2 binding or functional respiratory changes. This protocol defines experiments to test those unresolved hypotheses.
