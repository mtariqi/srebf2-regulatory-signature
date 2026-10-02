# Audit of 2026 EZH2–SREBP2 publication source data

## Source
Publication: An EZH2–SREBP2 axis promotes cholesterol biosynthesis
and represents a noncanonical vulnerability in tumorigenesis.
DOI: 10.1038/s41556-026-02048-x
Archive: 41556_2026_2048_MOESM12_ESM.zip
Audit date: 2026-10-02

## Inspection
The supplied archive contains 13 Excel workbooks and 10 PDFs.
Workbook sheets and cell contents were inspected.
PDF text was extracted to identify figure-panel contents.

## Findings
- Excel files contain figure measurements, statistical results,
  qPCR measurements, and selected gene-expression tables.
- PDFs contain source material for figure panels, including
  protein-blot panels.
- No BED, narrowPeak, broadPeak, or genomic peak-coordinate
  tables were identified.
- Some worksheets named gene_counts contain decimal expression
  values and selected gene subsets. Their normalization or
  transformation must be verified before reuse.
- The archive does not resolve the CUT&RUN control antibody
  identity or establish biological replication.

## Analysis decisions
- Do not use decimal expression tables as raw DESeq2 counts.
- Do not use selected gene subsets as a complete enrichment
  background.
- Do not treat this archive as a source of SREBP2 peak coordinates.
- Reserve these files for supporting context and figure comparison.
- Keep external expression results out of discovery signature
  selection when reserving that study for validation.
- Binding-supported signature construction remains pending
  suitable peak data or a justified analysis of raw binding reads.

## Next action
Audit GSE271000 for complete raw-count availability, direct
SREBP2 knockdown samples, matched controls, and biological
replication before selecting the validation analysis.
