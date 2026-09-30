# Research plan

## Candidate research gap
Evaluate reproducibility and transferability of an evidence-supported
SREBF2 transcriptional signature across independent datasets.
Confirm this gap through literature review before claiming novelty.

## Aims
1. Identify SREBF2-responsive genes using perturbation RNA-seq.
2. Integrate existing SREBP2 binding evidence.
3. Freeze the signature and validate it independently.

## Schedule
- Week 1: Literature review and dataset/sample audit.
- Week 2: Download inputs, record checksums, and perform QC.
- Week 3: Discovery differential-expression analysis.
- Week 4: Binding integration and signature construction.
- Week 5: Independent validation without refitting.
- Week 6: Compare with SREBF2 expression and cholesterol signatures.
- Week 7: Robustness and specificity analyses.
- Week 8: Figures, report, and reproducible release.
- Weeks 9–10: Buffer for raw-read processing or additional validation.

## Analysis rules
- Use raw integer counts for DESeq2, not TPM or FPKM.
- Verify biological replication and matched controls.
- Analyze studies separately and keep species separate.
- Select genes using discovery data only.
- Exclude SREBF2 itself from the target-based score.
- Freeze membership and scoring before validation.
- Report effect sizes, uncertainty, and adjusted significance.
- Treat binding plus perturbation support as candidate-target evidence.
- Do not infer protein activation or causality from RNA scores alone.
