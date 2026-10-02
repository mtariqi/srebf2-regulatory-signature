# Decision log

- Research analyses use real public data.
- Project focuses on SREBF2/SREBP2.
- Candidate datasets remain pending audit.
- Synthetic inputs are restricted to software tests.

Append dataset exclusions, contrast definitions, signature-selection
rules, and their supporting evidence here.

## Discovery cohort QC decision

- Retain all 12 samples with their original genotype and treatment labels.
- All samples support unstranded counting with featureCounts strand mode 0.
- SRR28966285 shows an atypical position relative to its KO-depletion replicates.
- Manifest and BAM read-group labels are consistent for SRR28966285.
- Label consistency does not independently authenticate biological identity.
- Primary analysis includes SRR28966285.
- Sensitivity analysis excludes SRR28966285 and compares effect estimates and signature stability.
- Retain SRR28966293; document elevated multimapping and consistent WT-FBS clustering.
- No sample is excluded or relabeled solely to improve clustering or significance.

## GSE267019 audit
- Human WT HeLa cells under cholesterol depletion.
- Assay measures H3K27ac, not SREBP2 occupancy.
- Two GEO samples listed: input GSM8258573 and H3K27ac GSM8258574.
- Genome assembly: hg38; listed processed files are BigWig tracks.
- Excluded from direct SREBP2 binding evidence.
- Retained as optional chromatin-context evidence, pending further audit.
- No replicated condition comparison established from these records.

## GSE282800 preliminary binding audit
- SREBP2 ChIP-seq in human WT HeLa cells; hg38 assembly.
- GSM8650746: FBS; GSM8650747: cholesterol depletion.
- One ChIP sample per condition listed; biological replication not established.
- No input-control sample listed within this series.
- Listed sample-level processed files are BigWig signal tracks.
- Supplementary archive contents and peak availability remain pending inspection.
- Candidate supporting occupancy evidence; not independent validation of discovery.
- No replicated differential-binding or direct-target claim established.

## GSE271001 preliminary audit
- SREBP2 CUT&RUN: HepG2 GSM8366948; Huh7 GSM8366956.
- Listed controls: HepG2 GSM8749017; Huh7 GSM8749018.
- Exact control type and suitability remain pending verification.
- One SREBP2 sample per cell line listed; biological replication not established.
- Processed assembly is hg19; cannot directly intersect with GRCh38 coordinates.
- Supplementary inventory confirms 19 BigWig files and no deposited peak-coordinate files.
- Overall design mentions 22Rv1, but no corresponding sample appears in the extracted records.
- Reserve as external binding-support candidate; do not use for signature selection yet.

## GSE324560 binding audit
- Human 22Rv1 and MV4;11 SREBP2 CUT&RUN; processed assembly hg38.
- Each cell line has Abcam and Cayman SREBP2 antibody samples and a listed control.
- Control antibody identity and biological replication remain unverified.
- Two antibodies are not automatically biological replicates.
- Supplementary inventory contains eight BigWigs and no peak-coordinate files.
- Reserve for external support; no signature selection from this study yet.
