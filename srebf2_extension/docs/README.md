# SREBF2 — Independent RNA-seq validation, STRING associations and SREBP2 occupancy

**Status (18 November 2026):** This is a reproducible *evidence-development package*, not three completed mechanistic validation experiments. It contains three rigorously labeled figures, scripts, assay/dataset inventories, and data templates. **Only primary HeLa gene-level effects have been calculated in Figure A.** Figure B intentionally contains no STRING edges until the API is successfully queried, and Figure C explicitly displays untested gene–assay pairs until peak overlaps are computed.

## Research questions and methods

1. **External RNA-seq validation:** Compare prespecified SREBF2 genetic perturbations (KO vs WT in HAP1 GSE147770; shSREBF2 vs shControl in HepG2 GSE271000) with original HeLa GSE267018 changes for the 12-gene cholesterol core. Recompute DESeq2 models from appropriate counts and replicates, select the correct treatment contrasts, harmonize gene identifiers, and compare directionality, effect size and multiplicity-adjusted significance. Statins and non-SREBF2 knockouts are *mechanistic comparators*, not independent SREBF2 genetic replications.
2. **Five-gene STRING association network:** Query human (taxonomy 9606) FDFT1, SQLE, CYP51A1, NSDHL, DHCR7 with no added interactors, required score 700/1000, functional network type, and record provenance. STRING associations combine predicted and curated evidence, are not necessarily direct physical binding, and cannot establish SREBP2 transcriptional occupancy. Only returned edges may be drawn.
3. **SREBP2 occupancy evidence:** Evaluate target promoters (SQLE, FDFT1, CYP51A1, DHCR24, LSS) in HeLa GSE282800 SREBP2 ChIP-seq and HepG2 ENCSR000EZO SREBF2 ChIP-seq; separately evaluate H3K27ac chromatin context in GSE267019. Require validated experimental peak files, controls, genome-assembly consistency, and strand-aware hg38 promoter coordinates (TSS -2,000 to +500 bp). ENCODE legacy HepG2 files may require hg19->hg38 liftOver, performed and audited before overlap. Sample conditions and cell lines differ; peak proximity cannot establish direct regulatory causality. Do not substitute H3K27ac peaks for SREBP2 occupancy.

## Verified source inventory

`data/dataset_screen.tsv` records 11 confirmed or index-reported datasets and how they may be used. GEO/ENCODE links are attached to each row. Further candidate data include CRISPR SREBF2 transcriptome studies reported in *Genome Biology* (2025), whose GEO series identifier requires accession-level verification before use. `GSE221582` is an **enzyme-KO**, not SREBF2-KO, microarray experiment; `GSE327534` is a **mouse Leydig-cell** model (secondary index).

## Outputs and evidence status

| Figure | Files | What is measured or illustrated | What is missing |
|---|---|---|---|
| A | `figures/Figure_A_External_Validation.{png,svg}` | Primary HeLa 12-gene log2FC from archived Phase15H. External cells display NA. | Gene-level effect estimates and padj from prespecified HepG2/HAP1 contrasts. |
| B | `figures/Figure_B_STRING_Network.{png,svg}` | Five named nodes to be queried; no fabricated associations. | True STRING API results and confidence-score edges. |
| C | `figures/Figure_C_Binding_Evidence_Audit.{png,svg}` | Assay availability and five high-priority target genes; status 'Not tested'. | Actual validated peaks, TSS intervals, controls, assembly harmonization and target overlaps. |

## Reproduce the three baseline figures

```bash
python3 -m pip install pandas numpy matplotlib networkx
python3 scripts/build_three_figures.py --data data --out figures
```

### Populate true STRING edges (requires network access)

```bash
python3 scripts/fetch_string_network.py --score 700 --out data/ppi_edges.tsv
python3 scripts/build_three_figures.py --data data --out figures
```

If the endpoint cannot be reached or returns no edges, the script exits without overwriting the existing edge table. Archive the raw API response and provenance JSON (automatically saved). Check identifier resolution and study bias before interpretation. Data are not available in this package because the STRING API was unreachable in the artifact-building environment.

### Populate independent RNA-seq data (requires independently computed contrasts)

Produce a TSV with the exact columns shown in `data/external_validation_template.tsv`, selecting **one prespecified contrast for each dataset** in Figure A. Positive log2FC must always mean increased in perturbed versus control cells. For GSE147770 use an explicitly identified +FBS and/or -FBS contrast. A HeLa sterol-depletion experiment and HAP1 serum starvation are NOT biologically equivalent.

```bash
cp data/external_validation_template.tsv data/external_validation.tsv
# Add real gene-wise DESeq2 contrasts (do not invent values).
python3 scripts/validate_external.py --external data/external_validation.tsv
python3 scripts/build_three_figures.py --data data --out figures
```

If analyzing both HAP1 contrasts, create separate figure runs or filter to one documented contrast; the plotting script fails if multiple contrasts for the same dataset are supplied simultaneously. Use original counts/appropriate design, not values transcribed from manuscript plots.

### Populate peak-level occupancy (requires harmonized real peak BED)

First obtain trustworthy peak files from each original dataset; check original peak-calling methods, actual genome build, controls and sample QC. Create `tss_hg38.tsv` with `gene`, `chrom`, `tss`, `strand` (1-based hg38 TSS), matching the original annotation definition, and obtain peak BED/narrowPeak files in hg38 (0-based coordinates). Then:

```bash
python3 scripts/check_peak_overlap.py \
  --tss tss_hg38.tsv --peaks /path/to/hela_SREBP2_hg38.narrowPeak.gz \
  --assay GSE282800_HeLa_SREBP2 --out data/hela_srebp2_overlap.tsv
python3 scripts/merge_occupancy_calls.py \
  --mapping GSE282800_HeLa_SREBP2=data/hela_srebp2_overlap.tsv \
  --out data/occupancy_targets_verified.tsv
# Back up original template before replacing:
cp data/occupancy_targets.tsv data/occupancy_targets_unverified.tsv
cp data/occupancy_targets_verified.tsv data/occupancy_targets.tsv
python3 scripts/build_three_figures.py --data data --out figures
```

`check_peak_overlap.py` is a coordinate-overlap audit, not a fully validated ChIP-seq analysis: peak reproducibility, specificity, input/background controls, statistical binding tests, chromosome naming and alternative TSSs still require scrutiny. Do not label missing peaks "not bound" without adequate assay sensitivity.

## Candidate captions (do not use figures A–C as completed validation results yet)

**A.** *SREBF2-dependent cholesterol-core expression: primary discovery and prespecified external validation framework.* Primary HeLa log2FCs are shown; NA indicates that the external dataset's corresponding gene-level test has not been supplied to this pipeline. No external replication is inferred.

**B.** *Five-gene sterol-biosynthesis STRING network: evidence-retrieval framework.* Nodes represent prespecified genes. Edges appear only when retrieved from STRING at the stated score cutoff; STRING associations are functional networks and cannot establish direct physical or transcriptional regulation.

**C.** *Public SREBP2 chromatin-occupancy audit for five cholesterol-biosynthesis genes.* Cells indicate the gene-level promoter-peak overlap validation status; the presence of a public ChIP-seq experiment does not establish occupancy at a specified gene. H3K27ac is a distinct chromatin mark.

## Manuscript limitation statement

"External SREBF2 perturbation experiments and public SREBF2 ChIP-seq assays were identified, but differences in cell type, sterol conditions, perturbation modalities, and assay controls preclude treating the studies as directly matched replications. Gene-level validation will require harmonized differential-expression contrasts, while binding claims will require peak-level analysis with experimentally appropriate controls and genome-coordinate harmonization. STRING associations represent functional evidence that does not necessarily imply direct physical interactions, and promoter proximity alone does not establish transcriptional regulation. Until those analyses have been completed and independently checked, these resources are presented as candidate validation datasets rather than confirmatory findings."

## Selected dataset and method references

- Xu et al. (2025), *Nature Metabolism*. [Publication](https://www.nature.com/articles/s42255-025-01291-0). Original GSE267018, GSE267019 and GSE282800.
- [GSE271000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271000): HepG2 SREBP2 shRNA RNA-seq series and other treatments.
- [GSE147770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147770): HAP1 SREBF2 KO vs WT with serum manipulation.
- [GSE113546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113546): four statins in human iPSC-cardiomyocytes.
- [GSE196187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196187): macrophage cholesterol reduction and inflammation conditions.
- [GSE221582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221582): HepG2 cholesterol-pathway enzyme knockout microarray.
- [GSE327534](https://www.omicsdi.org/dataset/geo/GSE327534): mouse Leydig Srebf2 ablation; secondary-index record.
- [ENCSR000EZO](https://www.encodeproject.org/experiments/ENCSR000EZO/): HepG2 SREBF2 ChIP-seq after pravastatin; [ENCSR000EZR](https://www.encodeproject.org/experiments/ENCSR000EZR/) control ChIP-seq.
- [STRING API documentation](https://string-db.org/help/api/) and [STRING methodology](https://doi.org/10.1093/nar/gkac1000).
- [Genome Biology 2025 HepG2 cholesterol screen and SREBF2 KO study](https://doi.org/10.1186/s13059-025-03531-8): verify accession before adding as an independent cohort.

**Do not overwrite the project's existing README or commit this package's placeholders as if they were measured biological results.** This folder is intended to be added as a new supplementary module, with figures included in the manuscript only after verification.
