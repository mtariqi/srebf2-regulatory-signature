# SREBF2 Regulatory Signature — Research Progress Report

**Reporting date:** 12 October 2026  
**Project:** Integrative Multi-Omics Investigation of SREBF2 Regulatory Responses  
**Repository:** https://github.com/mtariqi/srebf2-regulatory-signature  
**Researcher:** Md. Tariqul Islam  
**Collaborator:** Sarah Daley  
**Program:** MSc Bioinformatics, Northeastern University  
**Report type:** Cumulative research findings and 12 October documentation/repository progress  
**Status:** Transcriptomic analyses completed; exploratory epigenomic promoter-signal integration completed; proteomic/lipidomic integration proposed; computational manuscript in preparation.

> **Reporting convention.** The numerical findings below were generated during preceding research phases and are consolidated here as the evidence base available on 12 October 2026. The documented activity on this date is research synthesis, README revision, collaborator attribution, and GitHub version control. No new wet-lab assays, proteomics, lipidomics, or independent validation experiment is claimed for today. Values should be checked against the versioned result files before manuscript submission.

## 1. Executive summary

The project investigates how loss of **SREBF2**, encoding the sterol-responsive transcription factor **SREBP2**, affects transcriptional programs across cell lines and environmental contexts, and how those responses can be interpreted using complementary omics evidence. SREBF2 transcript abundance alone is not equivalent to nuclear, processed SREBP2 activity. A perturbation-derived expression signature is also not automatically a direct-target or binding-supported signature.

The HeLa discovery dataset **GSE267018** produced a frozen **1,326-gene** perturbation-derived candidate signature. In independent HepG2 shRNA data (**GSE271000**), the signature changed in the expected direction for **two of three** shRNA groups, whereas SREBF2 transcript abundance decreased in all three. In HAP1 SREBF2-knockout data (**GSE147770**), the aggregate signature response changed direction between serum-replete and serum-starved conditions. A factorial DESeq2 analysis identified **5,308** genes with a significant genotype-by-serum Wald interaction (BH FDR < 0.05). Of the mapped signature genes, **181** met both an expected directional-reversal criterion and interaction FDR < 0.05.

Enrichment of the 181-gene set highlighted **HALLMARK_OXIDATIVE_PHOSPHORYLATION** (13 overlapping genes; fold enrichment **5.246426**; BH-adjusted p = **6.149296 × 10⁻⁵**). HeLa SREBP2 ChIP-seq-derived promoter bigWig signals were detectable for all 13 genes, but expression- and biotype-matched control comparisons did **not** establish a consistent promoter-signal advantage. Consequently, the evidence supports **context-dependent SREBF2-associated transcriptional remodeling**, not demonstrated direct SREBP2 regulation of mitochondrial function.

The instructor's multi-omics requirement is addressed through **completed transcriptomics**, **exploratory epigenomic signal integration**, and explicit **public-data-first proteomics and lipidomics analysis plans**. Proposed methods are not represented as completed measurements.

## 2. Research rationale and questions

SREBP2 is a key regulator of sterol synthesis and uptake, but its activity is modulated by post-translational processing and cellular sterol status. Consequently, changes in SREBF2 mRNA, downstream transcription, protein abundance, and lipid phenotype may be discordant. Cell type, perturbation modality, and serum/sterol conditions can further modify the observed response.

**Primary research question:** To what extent are SREBF2-associated transcriptional responses transferable across cell types and perturbation conditions, and can public epigenomic evidence help prioritize plausible regulatory candidates?

**Secondary multi-omics question:** Can future proteomic and lipidomic evidence clarify whether the context-dependent transcriptomic signals propagate to protein-level or metabolic phenotypes?

**Working hypotheses:**

1. SREBF2 perturbation elicits both conserved and context-dependent gene-expression responses.
2. A frozen transcriptomic signature may capture some perturbation responses but may not consistently outperform SREBF2 RNA expression.
3. A subset of condition-reversal genes may be associated with mitochondrial/respiratory pathways, without necessarily being direct SREBP2 targets.
4. Orthogonal chromatin occupancy, protein, and lipid measurements could clarify mechanism; their incremental predictive value remains untested.

## 3. Data provenance and omics status

| Accession | Model / assay | Design or role | Evidence status |
|---|---|---|---|
| [GSE267018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267018) | HeLa RNA-seq | WT versus SREBF2 KO; FBS versus cholesterol depletion; 12 samples | Discovery and sensitivity analyses completed |
| [GSE271000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271000) | HepG2 RNA-seq | Controls and three shSREBP2 hairpins; eight samples | External scoring/diagnostic analyses completed |
| [GSE147770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147770) | HAP1 RNA-seq | WT versus SREBF2 KO; +FBS versus −FBS; 12 analyzed samples | External validation, factorial interaction, reversal audit and enrichment completed |
| [GSE282800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282800) | HeLa SREBP2 ChIP-seq-derived coverage | Promoter bigWig signal summary | Exploratory integration completed; no verified gene-specific peaks claimed |
| [GSE267019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267019) | HeLa H3K27ac ChIP-seq | Chromatin activity context | Dataset audited; integration not established as complete |
| [GSE271001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271001) | HepG2 SREBP2 CUT&RUN | Potential orthogonal occupancy resource | Candidate resource; no completed peak-based direct-target classification claimed |
| To be selected | Proteomics, LC-MS/MS | Differential protein abundance and RNA–protein concordance | Proposed; no analyzed cohort claimed |
| To be selected | Lipidomics, LC-MS | Sterols/lipid species and RNA–lipid pathway concordance | Proposed; no analyzed cohort claimed |

**Critical comparability restriction:** HAP1 −FBS reflects approximately four hours of serum starvation in the analyzed experiment. It is **not equivalent** to the HeLa cholesterol-depletion treatment. Cross-cohort comparisons are observational and do not isolate cell-line versus treatment effects.

## 4. Completed transcriptomic work

### 4.1 HeLa discovery: GSE267018

**Design:** Three biological samples per genotype-by-treatment group; 12 samples total. The DESeq2 model included genotype, treatment, and genotype-by-treatment interaction. Annotation used GENCODE v47; **17,158** genes were tested. Multiple testing used Benjamini–Hochberg correction.

| Contrast | FDR-significant genes | Increased in KO | Decreased in KO |
|---|---:|---:|---:|
| KO versus WT under cholesterol depletion | 6,387 | 3,105 | 3,282 |
| KO versus WT under FBS | 4,198 | 1,914 | 2,284 |
| Genotype × treatment interaction | 62 | 0 | 62 |
| Depletion contrast, sensitivity excluding `SRR28966285` | 6,835 | 3,196 | 3,639 |
| FBS contrast, same sensitivity | 4,671 | 2,147 | 2,524 |
| Interaction, same sensitivity | 175 | 28 | 147 |

These DEG counts use **BH FDR < 0.05 without an additional fold-change threshold**. The interaction count changed materially when `SRR28966285` was excluded; this must be retained as a sensitivity result, not substituted for the prespecified complete-sample analysis without justification.

### 4.2 Frozen discovery signature

The frozen `discovery_v1` signature was selected using the depletion KO-versus-WT contrast (adjusted p < 0.05, absolute log2 fold change ≥ 1), a concordant nonzero direction in the FBS and sample-exclusion sensitivity contrasts, exact annotation matching, and exclusion of SREBF2 itself. Statistical significance was **not** separately required in the FBS or sensitivity contrasts.

| Component | Genes | Direction under SREBF2 knockout |
|---|---:|---|
| Positive | 795 | Decreased expression |
| Negative | 531 | Increased expression |
| **Total** | **1,326** | Frozen candidate perturbation signature |

This is a **perturbation-derived**, not yet **binding-supported**, gene signature. The selection and scoring policies should remain frozen when reporting independent evaluations; validation-informed revisions must be labeled exploratory and tested in a new held-out cohort.

### 4.3 HepG2 independent validation: GSE271000

Eight samples comprised two controls and two replicates for each of three shRNA groups. All **1,326/1,326** signature gene IDs matched the HepG2 count matrix. A frozen direction-aware expression-percentile score was applied after count normalization.

| Group | Mean score | Difference from control | SREBF2 log2-normalized expression difference |
|---|---:|---:|---:|
| Control | 0.098951 | — | — |
| shSREBP2_664 | 0.089487 | −0.009463 | −1.680200 |
| shSREBP2_665 | 0.099560 | +0.000609 | −1.255098 |
| shSREBP2_667 | 0.088513 | −0.010438 | −2.474415 |

**Interpretation:** The signature decreased in two of three shRNA groups but not `shSREBP2_665`; SREBF2 RNA decreased in all three. With only two samples per group, these are primarily descriptive comparisons. The discordant hairpin must not be discarded post hoc. This dataset does **not** establish superior sensitivity, specificity, or accuracy of the signature relative to SREBF2 RNA. Library-layout discrepancies recorded in the validation audit should be retained in the methods/QC record.

### 4.4 HAP1 external validation: GSE147770

Twelve samples represent WT and SREBF2 KO under +FBS and −FBS, with three replicates per group. The discovery signature mapped to **1,065** genes (676 positive and 389 negative), rather than all 1,326 genes.

**Model:** DESeq2 `~ genotype + serum + genotype:serum`, with WT and +FBS as reference categories in the recorded analysis. **15,363** genes were tested.

| Analysis | BH FDR < 0.05 |
|---|---:|
| KO versus WT, +FBS | 3,954 |
| KO versus WT, −FBS | 7,729 |
| Genotype × serum interaction, Wald | **5,308** |
| Genotype × serum interaction, likelihood-ratio test | **5,309** |

Wald and likelihood-ratio interaction results had near-complete overlap. This is **methodological concordance on the same data**, not replication in another cohort. Across four scoring methods, the aggregate signature's KO-minus-WT response was negative under +FBS and positive under −FBS. The positive signature component exhibited the clearest reversal. This does not imply SREBP2 protein activity increased following knockout.

### 4.5 Directional-reversal audit

The audit combined condition-specific KO-versus-WT effect directions with the significance of the genotype-by-serum interaction.

| Component | Mapped | Both log2FC values evaluable | Expected reversal | Interaction FDR < 0.05 | Both reversal and interaction FDR < 0.05 |
|---|---:|---:|---:|---:|---:|
| Positive | 676 | 573 | 175 | 270 | **131** |
| Negative | 389 | 284 | 82 | 89 | **50** |
| **Combined** | **1,065** | **857** | **257** | **359** | **181** |

The 181 genes are a **selected, interaction-supported reversal subset**; they should not be presented as all SREBF2 targets or as genes directly bound by SREBP2.

## 5. Functional enrichment and mitochondrial candidates

### 5.1 Enrichment framework

Hallmark, Reactome, and GO Biological Process testing evaluated **5,647 pathways**, using the documented tested-gene background and BH adjustment. Results included:

- Positive 131-gene reversal subset: **16** FDR-significant pathways.
- Negative 50-gene reversal subset: **0** FDR-significant pathways.
- Combined 181-gene reversal subset: **7** FDR-significant pathways.
- Top combined pathway: **HALLMARK_OXIDATIVE_PHOSPHORYLATION**; overlap **13 genes**; fold enrichment **5.246426**; BH-adjusted p **6.149296 × 10⁻⁵**.

Additional significant pathways concerned respiratory electron transport, aerobic respiration, ATP synthesis, and mitochondrial translation. Enrichment remained evident under alternative tested backgrounds, although some sensitivity p-values were **nominal** rather than pathway-wide FDR-adjusted. No significant lipid/stress pathways in this selected 181-gene subset does **not** imply that cholesterol-related transcriptional changes were absent in the complete datasets.

### 5.2 Thirteen OXPHOS-associated reversal genes

`NDUFB7`, `ETFB`, `COX6B1`, `POR`, `MRPS12`, `MRPL34`, `NDUFA2`, `TIMM10`, `NDUFS6`, `ATP5MC3`, `ATP5ME`, `UQCRH`, `ATP5MF`.

All 13 belonged to the **positive** frozen signature component, exhibited negative KO-versus-WT log2 fold change in +FBS and positive KO-versus-WT log2 fold change in −FBS, and had a genotype-by-serum interaction at FDR < 0.05. The genes are **transcriptional candidates associated with an enriched pathway**, not proof of altered respiratory flux, ATP synthesis, mitochondrial function, or direct transcription-factor regulation.

### 5.3 Phase 6I candidate prioritization

The candidate ranking combined equal-weight within-set percentile ranks for RNA interaction significance (`−log10(padj)`) and reversal magnitude. **Promoter signal was not used in the ranking.**

| Rank | Gene | RNA-only relative priority score |
|---:|---|---:|
| 1 | **ETFB** | 100.00 |
| 2 | **TIMM10** | 88.46 |
| 3 | **MRPS12** | 76.92 |
| 4 | **NDUFB7** | 73.08 |
| 5 | **ATP5MC3** | 73.08 |

Scores rank these 13 candidates relative to one another; they are not posterior probabilities, functional validation scores, or measures of direct binding.

## 6. Epigenomic integration and evidence boundaries

### 6.1 Completed exploratory analysis

Public **GSE282800** HeLa SREBP2 ChIP-seq-derived bigWig signal was summarized over promoter windows approximately **TSS −2,000 bp to +500 bp**. For the 13 OXPHOS-associated candidates:

- **13/13** had detectable nonzero promoter coverage in at least one analyzed track.
- **10/13** had higher descriptive promoter signal in the cholesterol-depletion condition than FBS.
- A broad, unmatched signature-background comparison found higher median FBS promoter signal for the 13-gene set (exploratory unadjusted p approximately **0.011**).
- The unmatched promoter condition-difference comparison did not support enrichment (p approximately **0.658**).
- Crucially, **500 expression- and biotype-matched control draws did not show a consistent promoter-signal advantage**.

### 6.2 What cannot be concluded

A bigWig signal track is not a called, reproducible SREBP2 peak set. Nonzero promoter signal does not establish a gene-specific SREBP2 binding event, differential occupancy, or causal target regulation. Moreover, HeLa cholesterol-depletion promoter signal is not experimentally matched to HAP1 serum-starvation RNA effects. The promoter results should be presented as **exploratory cross-study evidence**, with matched-null findings emphasized rather than omitted.

### 6.3 Proposed epigenomic strengthening

If suitable public raw or processed occupancy data and controls are available: validate assay metadata and genome assembly; inspect replicate QC and input/negative controls; call or obtain reproducible SREBP2 peaks; use prespecified promoter/peak-to-gene assignment rules; perform matched-background enrichment and sensitivity analyses. **GSE267019 H3K27ac** may help contextualize active chromatin but is **not** SREBP2-specific occupancy evidence. **GSE271001** may provide additional CUT&RUN context if appropriate files and controls are available. These steps remain future work unless separately documented as completed.

## 7. Instructor-required multi-omics framework

The project must incorporate multiple omics **methods** while clearly distinguishing the layers for which quantitative data have actually been analyzed.

| Layer | Proposed or completed workflow | Principal output | Status |
|---|---|---|---|
| Transcriptomics | RNA-seq QC, normalization, DESeq2 factorial models, frozen signatures, enrichment | Gene effects, interaction FDR, candidate scores | **Completed** for analyzed cohorts |
| Epigenomics | SREBP2 ChIP-seq-derived promoter coverage; prospective peak-based occupancy integration | Promoter-signal audit; future reproducible peaks | **Exploratory promoter analysis completed**; peak-level classification proposed |
| Chromatin activity | H3K27ac ChIP-seq analysis, promoter/enhancer context | Differential chromatin-associated signals | **Proposed beyond dataset audit** |
| Proteomics | LC-MS/MS, peptide/protein QC, normalization, differential abundance, pathway mapping | Protein-level response and RNA–protein concordance | **Proposed** |
| Lipidomics | LC-MS lipid/sterol profiling, internal standards, pooled QC, differential lipid abundance | Lipid-level response and pathway concordance | **Proposed** |
| Integrated modeling | Cross-omics feature mapping, evidence tiers, benchmark against simple baselines | Independently assessed generalizability | **Proposed** |

### 7.1 Public-data-first proteomics plan

Search repositories such as PRIDE/ProteomeXchange for **SREBF2/SREBP2 perturbation** proteomic datasets with documented cell type, perturbation, treatment, replicate structure, and accessible protein quantification. If a compatible dataset exists: perform contaminant filtering and QC; assess missingness and batch structure; normalize appropriately; fit a prespecified protein-level model; adjust for multiple testing; harmonize protein identifiers to gene symbols; test directional concordance with RNA and the 13 candidates. Cross-study associations should be reported as descriptive unless experiments are matched. No protein abundance result is reported at present.

### 7.2 Public-data-first lipidomics plan

Search Metabolomics Workbench, MetaboLights, and relevant supplementary datasets for compatible sterol/lipid profiles. Record sample provenance, treatment, internal standards, lipid annotations, instrument/batch metadata, and QC samples. Where data quality permits, assess total/free cholesterol, cholesterol esters, and related lipid classes; perform suitable normalization and multiplicity correction; integrate pathway-level trends with RNA/protein evidence. No lipid phenotype is claimed as measured or validated in this repository at present.

### 7.3 Cross-omics evidence hierarchy

Use an explicit evidence matrix: (i) RNA interaction and reversal; (ii) promoter signal only; (iii) reproducible occupancy peak, if obtained; (iv) independent protein response, if obtained; (v) lipid/functional phenotype, if obtained. Missing layers should remain missing, not be imputed as positive support. Compare any new integrated activity score against the unchanged 1,326-gene score, SREBF2 RNA expression, and prespecified sterol-pathway scores in **new held-out data**. Because existing validation cohorts have been inspected, a redesigned model tested on them alone would be exploratory.

## 8. Proposed Phase 6J validation (not performed)

The documented future protocol uses a genotype × serum factorial design in HAP1 or another suitable model, with WT and SREBF2 KO under +FBS and −FBS. Priorities include confirming treatment duration and biological replicate definitions; RT-qPCR for `ETFB`, `TIMM10`, `MRPS12`, `NDUFB7`, and `ATP5MC3`; a prespecified interaction analysis with multiple-testing control; and, if resources become available, orthogonal SREBP2 occupancy, rescue, and mitochondrial respiration measurements. The protocol is a future experimental roadmap and does **not** represent performed assays. See `protocols/phase6j/README.md`.

## 9. Computational methodology and reproducibility

**Software/methods recorded in the project:** R/DESeq2 for count-based negative-binomial models; Python for scoring, identifier mapping, promoter-signal processing, and candidate ranking; BH false-discovery control; GENCODE v47 for discovery annotation; versioned output tables and manifests. Environment details recorded during the project include R 4.6.1, DESeq2 1.52.0, and Bioconductor 3.23. Verify exact runtime versions against manifests before publication.

**Representative scripts and artifacts** (paths should be checked against the current repository before claiming a clean fresh-clone reproduction):

| Purpose | Representative repository path |
|---|---|
| Discovery summary | `results/de/discovery_v1/DE_SUMMARY.tsv` |
| Frozen signature | `results/signature/discovery_v1/` |
| HepG2 scores | `results/signature_validation/GSE271000_v1/` |
| HAP1 interaction | `results/validation/GSE147770/interaction_v1/` |
| HAP1 reversal audit | `results/validation/GSE147770/direction_audit_v1/` |
| Enrichment | `results/validation/GSE147770/enrichment_v1/` |
| RNA candidate ranking | `results/validation/GSE147770/candidate_prioritization_v1/` |
| HAP1 interaction script | `scripts/srebf2_GSE147770_phase6d_interaction_v1.R` |
| Reversal audit script | `scripts/srebf2_GSE147770_phase6e_direction_audit_v1.R` |
| Enrichment script | `scripts/srebf2_GSE147770_phase6f_enrichment_v1.R` |
| Promoter audit script | `scripts/srebf2_GSE147770_phase6h_binding_audit_v1.py` |
| Candidate ranking script | `scripts/srebf2_GSE147770_phase6i_candidate_prioritization_v1.py` |
| Proposed validation protocol | `protocols/phase6j/README.md` |

**Reproducibility safeguards:** Retain frozen signature and scoring policies; store hashes/manifests; separate raw sequencing reads and large intermediates from Git; distinguish exploratory analyses from frozen evaluation; report all hairpins and sample-exclusion sensitivities; record annotation and genome-assembly compatibility; avoid claiming independent replication when tests reuse the same samples.

**Read-only checks (fish-compatible):**

```bash
cd ~/srebf2-regulatory-signature
cat results/de/discovery_v1/DE_SUMMARY.tsv
cat results/validation/GSE147770/interaction_v1/RESULT_SUMMARY_v1.tsv
cat results/validation/GSE147770/direction_audit_v1/DIRECTION_SUMMARY_v1.tsv
cat results/validation/GSE147770/enrichment_v1/ENRICHMENT_SUMMARY_v1.tsv
cat results/validation/GSE147770/candidate_prioritization_v1/EVIDENCE_SUMMARY_v1.tsv
```

These commands **inspect** existing results; they do not constitute a complete rerun. A publication-ready reproducibility release should specify input acquisition, sample manifests, checksums, environments, script order, and expected outputs.

## 10. Progress documented on 12 October 2026

### 10.1 README and scientific framing

A substantially expanded README was prepared to foreground the completed HAP1 validation, interaction/reversal findings, OXPHOS enrichment, exploratory promoter audit, candidate ranking, and the instructor-required four-omics research framework. The project is positioned as a **public-data-driven integrative multi-omics investigation**, not a completed four-assay experimental study. Manuscript preparation is the next major writing milestone.

### 10.2 Collaboration and citation

**Sarah Daley** was added as collaborator in the revised README and its suggested repository citation:

> Islam, M. T., & Daley, S. (2026). *SREBF2 regulatory signature: Integrative transcriptomic and epigenomic analysis with a proposed multi-omics framework* [Research code and data repository]. GitHub. https://github.com/mtariqi/srebf2-regulatory-signature

This is a **suggested software/repository citation**, not a published journal article. Manuscript authorship order and contribution statements require mutual agreement before submission.

### 10.3 GitHub activity and outstanding README publishing step

The user-provided terminal output confirms successful commit **`9cca0c4`**, message **`Add Sarah Daley as collaborator and update repository citation`**, pushed to `main`. That commit **created `README_SREBF2_UPDATED_V2.md`** (398 insertions). It did **not** itself replace the root `README.md`, which GitHub normally displays on the repository homepage. A later command sequence was recommended to copy the updated document into `README.md`; its successful execution has **not been verified** in the available terminal output. Therefore, this report does not claim that the root README has already been updated.

Earlier research results were pushed in commit `a3bd670` as part of the validation-results publication milestone. Separate GitHub Actions build checks were reported failing in a prior screenshot; the root cause has **not been diagnosed from job logs**. CI status is a software/reproducibility issue, not by itself evidence against the biological findings.

### 10.4 Working-tree hygiene and provenance

Three tracked validation result files were previously reported modified, potentially after whitespace cleanup:

- `results/validation/GSE147770/enrichment_v1/RUN_MANIFEST_v1.txt`
- `results/validation/GSE271000/background_v2/BACKGROUND_COMPARISON_v2.tsv`
- `results/validation/GSE271000/processed_DE_v1/VALIDATION_GENE_MAPPING.tsv`

Untracked promoter-analysis scripts, `results/binding/`, QC outputs, metadata, and an earlier progress report were also present. **Do not use `git add .` indiscriminately.** Inspect differences before committing or restoring; trailing TSV delimiters may encode empty fields. These working-tree observations describe the last shared terminal status and may have changed locally.

## 11. Interpretation and contribution

The strongest current scientific contribution is the observation that **SREBF2 perturbation-associated transcriptional responses are condition dependent**, including a substantial genotype-by-serum interaction and a direction-reversing subset enriched for mitochondrial oxidative phosphorylation-related genes. The frozen signature's variable transferability across HepG2 hairpins and HAP1 serum conditions is a meaningful negative/qualified result, not a reason to omit the signature analysis. The lack of robust matched-null promoter-signal enrichment prevents elevating the 13 candidates to verified direct SREBP2 targets.

The results support a cautious manuscript framing such as:

> **Context-Dependent SREBF2 Perturbation Responses Reveal a Mitochondrial Transcriptional Signature: An Integrative RNA-seq and Promoter-Signal Analysis**

This is a proposed working title, not a submitted or accepted publication.

## 12. Limitations and threats to inference

1. **Different biological systems:** HeLa, HepG2, and HAP1 differ in cell identity and experimental background.
2. **Different perturbations:** Knockout and shRNA knockdown cannot be assumed equivalent; shRNA off-target effects are possible.
3. **Non-equivalent treatment conditions:** HAP1 serum starvation differs from HeLa cholesterol depletion.
4. **Small validation groups:** HepG2 has two replicates per group, limiting inference.
5. **Signature transfer:** Only 1,065/1,326 genes mapped in HAP1; mapped composition differs from discovery.
6. **Interaction model sensitivity:** HeLa interaction discoveries changed following sample exclusion.
7. **Selection and multiplicity:** Reversal candidates were selected using both direction and significance; enrichment and exploratory comparisons require explicit backgrounds and adjusted p-values.
8. **Binding uncertainty:** Promoter bigWig coverage is not validated occupancy; matched-null analysis did not show consistent advantage.
9. **No functional validation:** No direct measurements of respiration, mitochondrial function, protein abundance, or sterol/lipid concentrations were completed as part of the reported analyses.
10. **No proven model superiority:** Neither a frozen signature nor a future multi-omics model has demonstrated consistent superiority to SREBF2 RNA expression across independent conditions.
11. **Public-data and reproducibility dependencies:** Dataset metadata, licenses, software versions, script availability, and CI failures require a final audit before publication.

## 13. Prioritized next steps

| Priority | Task | Expected deliverable | Completion criterion |
|---|---|---|---|
| P0 | Publish revised content as root `README.md` | Updated repository homepage | Verify commit and rendered README |
| P0 | Audit changed TSV/manifests and untracked analysis scripts | Clean, attributable working tree | Review diffs; selectively commit code/results |
| P0 | Inspect failing GitHub Actions logs | Minimal CI repair | Appropriate syntax/smoke checks pass |
| P1 | Freeze a manuscript analysis inventory | Source-to-figure/table map | Each reported number linked to a versioned output |
| P1 | Create manuscript outline, methods, results, discussion | Initial draft | Claims match documented evidence boundaries |
| P1 | Prepare publication figures and captions | Figure set | Readable labels, reproducible scripts, stated tests |
| P2 | Audit suitable public proteomics datasets | Feasibility matrix | Explicit inclusion/exclusion and assay compatibility |
| P2 | Audit suitable public lipidomics datasets | Feasibility matrix | Provenance, QC and metadata assessed |
| P2 | Assess peak-level SREBP2 occupancy data | Exploratory binding evidence tier | Controls/replicates/assembly documented |
| Future | Orthogonal Phase 6J experiments | RT-qPCR, occupancy and functional assays | Only if resources and collaboration permit |

**Near-term decision:** Prioritize reproducibility, instructor-facing multi-omics documentation, and manuscript drafting. New wet-lab experiments are not required to honestly report the completed computational results; their absence must be disclosed.

## 14. GitHub submission instructions

Save this report under `reports/SREBF2_Progress_Report_2026-10-12.md` (or another existing documentation folder). From the repository root:

```bash
mkdir -p reports
cp ~/Downloads/SREBF2_Progress_Report_2026-10-12.md reports/
git status --short
git add reports/SREBF2_Progress_Report_2026-10-12.md
git diff --cached --check
git diff --cached --stat
git commit -m "Document SREBF2 multi-omics research progress for 2026-10-12"
git push origin main
```

The `cp` example assumes the report was downloaded to `~/Downloads`; adjust the source path if necessary. The selective `git add` avoids staging unrelated modified or untracked research files.

## 15. Attribution and citation

**Researcher:** Md. Tariqul Islam  
**Collaborator:** Sarah Daley  
**Repository:** https://github.com/mtariqi/srebf2-regulatory-signature

**Suggested repository citation:**

Islam, M. T., & Daley, S. (2026). *SREBF2 regulatory signature: Integrative transcriptomic and epigenomic analysis with a proposed multi-omics framework* [Research code and data repository]. GitHub. https://github.com/mtariqi/srebf2-regulatory-signature

**Data attribution:** Cite the original GEO studies and their associated publications for GSE267018, GSE271000, GSE147770, GSE282800, and any subsequently used public datasets. Cite relevant software, annotation releases, pathway collections, and data repositories in the manuscript. The suggested repository citation is not a substitute for primary dataset attribution.

---

**Report conclusion:** The project has progressed from a HeLa perturbation signature to a cross-context transcriptomic investigation with an exploratory epigenomic evidence layer and a clearly specified multi-omics expansion plan. The key reproducible result is context-dependent SREBF2-associated transcriptional reversal with enrichment of oxidative phosphorylation-related genes; direct SREBP2 binding and downstream protein/lipid effects remain unverified. The immediate focus is a transparent, publication-oriented computational manuscript and repository reproducibility audit.
