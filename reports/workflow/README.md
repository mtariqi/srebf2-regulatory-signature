# Phase 12 Snakemake Reproducibility

This directory documents the initial Snakemake orchestration
and reproducibility testing of the SREBF2 Phase 12 analyses.

## Verification

- Snakemake version: 7.32.4
- Python version: 3.14.4
- Eight Phase 12 TSV files were checked.
- Six files matched their original SHA-256 checksums.
- Two files differed only by floating-point rounding:
  - PHASE12_EXPRESSION_STRATIFIED_MODELS.tsv
  - PHASE12_VIF.tsv
- Both mismatched files passed numerical equivalence testing
  with rtol=1e-8 and atol=1e-10.
- Maximum observed absolute difference: approximately 8.88e-16.

## Scope and limitations

This validation was conducted in the existing local software
environment using the previously generated regulatory-signal
matrix as input.

It does not establish reproducibility of upstream raw-data
processing or execution across different operating systems,
containers, or HPC environments.

## Files

- phase12_before.sha256: Original output checksums.
- phase12_numerical_reproducibility.tsv: Numerical comparison.
- phase12_reproducibility_summary.txt: Verification summary.
- phase12_reproducibility_verification.txt: Initial checksum test.
- phase12_dag.svg: Snakemake dependency graph.
- phase12_rulegraph.svg: Workflow rule graph.
