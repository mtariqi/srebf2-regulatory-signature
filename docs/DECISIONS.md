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
