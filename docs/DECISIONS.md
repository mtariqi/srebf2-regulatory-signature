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
