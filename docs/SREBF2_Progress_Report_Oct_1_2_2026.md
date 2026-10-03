# SREBF2 regulatory-signature project: combined progress report

**Prepared for:** project groupmate
**Reporting period:** October 1–2, 2026 (Toronto time)
**Status cutoff:** October 2, approximately 10:13 p.m.
**Repository:** https://github.com/mtariqi/srebf2-regulatory-signature

## Executive overview

I have completed the discovery RNA-seq processing and differential-expression analysis for 12 HeLa samples, annotated the results, and begun preparing an external HepG2 validation cohort. Both external controls now have verified gene-count tables. The first SREBP2 knockdown sample has passed download, archive, FASTQ-pair and compression checks; its FastQC flags have been reviewed and the alignment command has been supplied. Its alignment completion has not yet been reported.

This is a computational reanalysis of public data, not a new laboratory experiment. The intended outcome is a reproducible SREBF2 perturbation-response signature and an evaluation of its transfer to an external dataset. A final signature has not yet been frozen or externally validated. RNA changes alone do not establish SREBP2 protein activity or direct DNA binding.

This report reconstructs completed coding steps from the recorded commands and outputs. Some reference preparation began before October 1. Exact day attribution is not available for every discovery step, so those are reported as work completed by this reporting period rather than assigned an invented date. Code examples below summarize the actual workflow; placeholders are explanatory and are not a new executable pipeline.

## 1. Experimental design

### Discovery: GSE267018, HeLa

| Group | Samples | Purpose |
|---|---:|---|
| Wild type, FBS | 3 | Reference genotype under FBS |
| SREBP2 knockout, FBS | 3 | Genotype effect under FBS |
| Wild type, cholesterol depletion | 3 | Wild-type response under depletion |
| SREBP2 knockout, cholesterol depletion | 3 | Genotype effect under depletion |

The primary analysis includes all 12 samples. A predefined sensitivity fit excludes SRR28966285; it supplements the primary analysis rather than replacing it because it gives more significant genes.

### External expression cohort: GSE271000, HepG2

| Sample | Run | Label | Current status |
|---|---|---|---|
| GSM8366913 | SRR29927661 | shControl 1 | Alignment and counts validated |
| GSM8366914 | SRR29927662 | shControl 2 | Alignment and counts validated |
| GSM8366917 | SRR29634180 | shSREBP2 664, label 1 | FASTQ QC reviewed; alignment result pending |
| GSM8366918 | SRR29634179 | shSREBP2 664, label 2 | Processing pending |
| GSM8366919 | SRR29634178 | shSREBP2 665, label 1 | Processing pending |
| GSM8366920 | SRR29634177 | shSREBP2 665, label 2 | Processing pending |
| GSM8366921 | SRR29634176 | shSREBP2 667, label 1 | Processing pending |
| GSM8366922 | SRR29634175 | shSREBP2 667, label 2 | Processing pending |

The labels indicate two samples per hairpin, but biological independence must be confirmed from study documentation. The three hairpin comparisons share controls and are not three independent studies.

## 2. Reference preparation and integrity

I downloaded the GRCh38 primary-assembly genome and GENCODE v47 primary-assembly GTF, tested compressed-file integrity and recorded SHA-256 checksums.

```fish
curl --fail --location --retry 3 --output FILE URL
pigz -t data/raw/reference/gencode_v47/GRCh38.primary_assembly.genome.fa.gz
pigz -t data/raw/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf.gz
sha256sum FILES > metadata/gencode_v47_reference.sha256
```

A genome–annotation compatibility check examined 4,117,647 annotation features against 194 genome contigs and reported zero errors. Splice-site extraction exited successfully and produced 528,703 entries.

I built a genome-only HISAT2 index using one thread to manage memory usage. Annotated splice sites are supplied separately during alignment.

```fish
hisat2-build -p 1 GENOME_FASTA INDEX_PREFIX
hisat2-inspect -s INDEX_PREFIX
sha256sum INDEX_FILES > metadata/GRCh38_gencode_v47_hisat2_index.sha256
sha256sum SPLICE_SITE_FILE > metadata/gencode_v47_splicesites.sha256
```

The build took approximately 1 hour 3 minutes. Build and inspection exit statuses were zero. Validation confirmed eight nonempty index files and 194 reference sequences. Reference summaries and checksums were committed.

## 3. Discovery raw-read processing

For each discovery run, I downloaded the SRA archive, validated it, extracted FASTQs, compressed them, checked gzip integrity, recorded checksums, ran FastQC and checked mate correspondence.

```fish
prefetch RUN -O SRA_DIRECTORY --max-size SIZE
vdb-validate ARCHIVE
fasterq-dump ARCHIVE --split-files --threads 4 --outdir FASTQ_DIRECTORY --temp TEMP_DIRECTORY
pigz -p 4 MATE1.fastq MATE2.fastq
pigz -t MATE1.fastq.gz MATE2.fastq.gz
sha256sum MATE1.fastq.gz MATE2.fastq.gz > CHECKSUM_FILE
fastqc --threads 2 --outdir QC_DIRECTORY MATE1.fastq.gz MATE2.fastq.gz
```

The Python pair validator streamed both mates, checked four-line FASTQ structure, sequence/quality length equality, matching read identifiers, synchronized end-of-file and the expected pair total. This avoids loading the files into memory.

The remaining discovery samples were processed sequentially using `scripts/process_remaining_srebf2.py`. The queue completed with exit status zero; individual sample QC remained subject to review.

## 4. Alignment and counting

Reads were aligned using HISAT2 with four threads, the verified genome index, known splice sites and run/sample read groups. Output was streamed directly into samtools sort to avoid a large intermediate SAM file.

```fish
hisat2 -p 4 -x INDEX_PREFIX --known-splicesite-infile SPLICE_SITES \
    --rg-id RUN --rg SM:SAMPLE --rg PL:ILLUMINA \
    -1 MATE1.fastq.gz -2 MATE2.fastq.gz 2> HISAT2_LOG \
| samtools sort -@ 1 -m 768M -T TEMP_PREFIX -o OUTPUT.sorted.bam - 2> SORT_LOG
set codes $pipestatus
samtools quickcheck -v OUTPUT.sorted.bam
samtools index -@ 2 OUTPUT.sorted.bam
samtools flagstat -@ 2 OUTPUT.sorted.bam > FLAGSTAT.txt
```

Both pipeline exit statuses were checked. BAM integrity, indexing and flagstat supplied additional checks.

Gene counting was performed with featureCounts in each of three strand modes:

```fish
featureCounts -T 4 -p --countReadPairs -B -C -s MODE \
    -t exon -g gene_id -a GENCODE_GTF -o COUNT_TABLE SORTED_BAM
```

Counts aggregate exons by gene ID, count paired fragments, require both mates to be mapped and exclude chimeric pairs under these options. The logs indicated multimapping and multi-overlapping reads were not counted.

Python validation checked unique gene IDs, nonnegative integer counts and equality between the table sum and the summary Assigned count. Compact two-column gene-count tables were retained in Git; large full featureCounts annotation tables remain local. Removing the large tables in a later commit does not remove their earlier Git history.

## 5. Discovery cohort QC

All 12 discovery samples aligned at 97.49–98.09%. Unstranded assignment ranged from 67.88–75.48%; forward and reverse assignments were similar and lower, supporting unstranded counting.

SRR28966293 had lower assignment (67.88%) and higher multimapping (12.57%). This flag was recorded for review rather than used alone to exclude the sample. Its reported expression correlations with the other WT_FBS samples were high, but correlation alone cannot establish sample quality or identity.

FastQC repeatedly flagged sequence content and duplication, with read-1 tile warnings. These flags were reviewed alongside other evidence; they were not automatically converted into trimming, deduplication or exclusion decisions.

The merged count matrix passed validation: **78,932 genes × 12 samples**, all nonnegative integers. Assigned library totals ranged from 19,318,488 to 23,852,695. DESeq2 size factors ranged approximately 0.921–1.090.

`scripts/srebf2_sample_qc.R` generated library-size QC, PCA coordinates/figures, sample correlations and methods/session records. The image files were not successfully viewed in this conversation; no visual PCA conclusion is asserted here.

## 6. Analysis environment and differential expression

An initial R package-loading error arose from incompatible personal-library packages. I moved analysis into the isolated `srebf2-analysis` environment and confirmed R 4.5.3 and DESeq2 1.50.2 loaded successfully.

```fish
conda activate srebf2-analysis
Rscript --vanilla -e '
.libPaths(file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library"), include.site=FALSE)
stopifnot(requireNamespace("DESeq2", quietly=TRUE))
cat(R.version.string, "\n")
cat(as.character(packageVersion("DESeq2")), "\n")
'
```

Environment specifications and R session information were saved. `scripts/srebf2_differential_expression.R` fitted a genotype-by-treatment design with WT and FBS as reference levels. The fixed filter retained genes with at least 10 counts in at least three of the original 12 samples, yielding 17,158 genes for both fits.

The three comparisons were knockout versus wild type under depletion, knockout versus wild type under FBS, and their difference (interaction). The interaction is a difference of genotype effects, not simply a comparison of depletion samples.

| Fit | Contrast | Genes at adjusted p < .05 | Up | Down |
|---|---|---:|---:|---:|
| All 12 | KO vs WT, depletion | 6,387 | 3,105 | 3,282 |
| All 12 | KO vs WT, FBS | 4,198 | 1,914 | 2,284 |
| All 12 | Interaction | 62 | 0 | 62 |
| Excluding SRR28966285 | KO vs WT, depletion | 6,835 | 3,196 | 3,639 |
| Excluding SRR28966285 | KO vs WT, FBS | 4,671 | 2,147 | 2,524 |
| Excluding SRR28966285 | Interaction | 175 | 28 | 147 |

Both fits completed successfully. Some interaction adjusted p-values are unavailable; this must be distinguished from a nonsignificant adjusted p-value when comparing fits. Counts of significant genes alone do not show robustness, direct regulation or novelty.

## 7. Gene annotation and preliminary biological observations

`scripts/annotate_srebf2_results.py` mapped the three result tables from each fit and the three sensitivity tables to GENCODE v47. Every table had 17,158 rows, all mapped and none unmapped. Outputs include `ANNOTATION_QC.tsv` and `SREBF2_RESULTS.json`.

SREBF2 itself showed log2 fold changes of −2.719 under depletion and −2.723 under FBS in the primary fit, with very small adjusted p-values. These correspond to approximately 85% lower RNA abundance. They do not prove complete absence of protein. Its interaction was close to zero and not significant.

The supplied cholesterol-related gene panel showed consistent decreases under depletion. Examples include HMGCR (−3.180), SQLE (−3.454), MVD (−4.081) and LDLR (−1.526). Several had stronger genotype effects under depletion than FBS. Seven panel genes had significant interactions in both fits: DHCR24, FDPS, HMGCS1, HMGCR, SQLE, DHCR7 and MVD. FDFT1 and MVK crossed the significance threshold differently between fits. These are preliminary discovery observations, not a final selected signature or proof of direct binding.

## 8. Binding-data and external-data audits

I inspected GEO metadata and supplementary inventories before treating datasets as binding evidence.

| Dataset | Audit outcome | Intended use/limitation |
|---|---|---|
| GSE267019 | HeLa H3K27ac, not SREBP2 immunoprecipitation | Chromatin context; does not establish SREBP2 occupancy |
| GSE282800 | HeLa SREBP2 ChIP-seq, hg38; two deposited BigWigs | Supporting signal; no deposited peak files, replication/input not established within series |
| GSE271001 | HepG2/Huh7 SREBP2 CUT&RUN, reported hg19; BigWigs | External support candidate; assembly and control suitability need resolution |
| GSE324560 | 22Rv1/MV4;11 SREBP2 CUT&RUN, hg38; eight BigWigs | Two antibodies are not automatically biological replicates; controls unverified |
| GSE287355 | Superseries linking six component studies | Study inventory and provenance |

Source metadata were downloaded with curl, gzip-tested, hashed and relevant fields extracted with `rg` or Python. Findings and limitations were recorded in `docs/DECISIONS.md` and audit files. Source-data inspection did not establish a usable peak-coordinate resource or resolve all control details. BigWig signal is not automatically a validated peak set.

Three deposited GSE271000 hairpin tables were downloaded and inspected: 664 (17,682 rows), 665 (17,249), 667 (17,197). They contain DE statistics, not sample-level expression columns. Their baseMean values cannot reconstruct individual sample counts. They are reserved for later response-concordance analysis after signature freezing; they cannot currently supply sample-level scoring validation.

## 9. External-control processing completed today

The correct GSE271000 run table contained 36 runs; selection found all eight required experiments. An initially misidentified 30-run file belonged to the discovery study and was replaced. The selected eight archives total approximately 23.33 GB in reported metadata.

Both controls were labelled SINGLE in the run table, but archive inspection and full FASTQ pair validation demonstrated paired-read files. Original metadata were preserved and observed structure recorded separately.

| Metric | SRR29927661 | SRR29927662 |
|---|---:|---:|
| Validated pairs | 52,137,798 | 59,864,779 |
| Read length | 101 bases | 101 bases |
| Overall alignment | 95.17% | 95.57% |
| Properly paired | 91.07% | 91.55% |
| Unstranded assignment | 67.02% | 67.23% |
| Forward assignment | 40.15% | 40.24% |
| Reverse assignment | 40.02% | 40.15% |
| Assigned fragments | 38,970,905 | 44,885,109 |
| Nonzero genes | 29,955 | 30,550 |
| Validated gene rows | 78,932 | 78,932 |
| Alignment/sort elapsed time | 36m 51s | 42m 37s |

Archive extraction used SRA Toolkit 3.2.1 in `srebf2-sra-test`; alignment/QC used `srebf2-rnaseq`. Controls passed gzip integrity, full mate validation, BAM quickcheck/indexing and count validation. Both controls had tile/content/duplication flags but passed per-base quality and adapter modules. Similar flags do not establish their cause. Zero duplicates in flagstat means no records were marked as duplicates, not absence of PCR duplication.

A fish filename mistake produced literal braces around SRR29927662 in count filenames. The files were renamed without overwriting existing targets; counting did not need repeating. Future paths use separate variable/string segments instead of braces inside a quoted string.

## 10. First knockdown: progress tonight

SRR29634180 (GSM8366917, hairpin 664 label 1) downloaded successfully. `vdb-validate` returned zero; the five inspected spots each had two biological reads of 101 bases. Extraction returned zero and produced 41,852,139 spots / 83,704,278 reads.

Both FASTQs were compressed with pigz, gzip-tested and validated over **all 41,852,139 pairs**. SHA-256 checksums were recorded and FastQC completed successfully.

Both mates passed per-base quality, sequence quality, GC, N, length and adapter checks. They failed tile quality, base content and duplication and warned for overrepresented sequences. The listed overrepresented sequences were poly-A/poly-T, individually 0.132–0.284% in mate 1 and 0.177–0.186% in mate 2. These are the reported FastQC sequences, not an estimate of every low-complexity read.

Mate 1 had 32/701 tiles with any deviation below −5; mate 2 had 72/701. These fractions describe tiles, not the fraction of affected reads. QC details were saved. Alignment was proposed with the same reference and settings as the controls; completion remains unconfirmed at this cutoff.

## 11. Coding and output inventory

| Code or recorded operation | Function | Output location |
|---|---|---|
| curl/gzip/pigz/sha256sum | Download, integrity and provenance | `metadata/`, `metadata/source/` |
| Compatibility checker | Genome/GTF contig and coordinate checks | `results/qc/reference/` |
| hisat2-build / hisat2-inspect | Reference index creation and inspection | `data/raw/reference/gencode_v47/`, QC summaries |
| SRA Toolkit prefetch/vdb-validate/vdb-dump/fasterq-dump | Archive acquisition, structure and extraction | `data/raw/`, `logs/` |
| Python streaming FASTQ validator | Record and mate validation | `PAIR_VALIDATION.txt` |
| FastQC plus Python ZIP parser | Read QC and flag extraction | Per-sample `QC_SUMMARY.txt`, `QC_FLAG_DETAILS.txt` |
| HISAT2/samtools | Alignment, sorting, BAM integrity/index/QC | `data/processed/`, alignment summaries |
| `scripts/count_pilot_strand_check.fish` | Discovery pilot strand comparison | Pilot strand summaries |
| `scripts/process_remaining_srebf2.py` | Sequential discovery processing | Counts and per-run QC |
| Python count validator/compact exporter | Integer/ID/total checks | `results/counts/`, `COUNT_VALIDATION.txt` |
| Python count-matrix assembly | Merge validated per-run counts | `GSE267018_raw_count_matrix.tsv` |
| `scripts/srebf2_sample_qc.R` | Expression-level QC | `results/qc/cohort/expression/` |
| `scripts/srebf2_differential_expression.R` | Two DESeq2 fits and three contrasts | `results/de/discovery_v1/` |
| `scripts/annotate_srebf2_results.py` | GENCODE result annotation | `results/de/discovery_v1/annotated/` |
| Python/rg GEO audits | Sample, assay, assembly and file checks | `metadata/` audits and decisions |
| git add/diff/commit/push | Versioned checkpoints | GitHub repository |

Some validation operations are currently inline Python commands or interactive fish functions, not confirmed standalone tracked scripts. Consolidating them into a saved, resumable validation workflow is an outstanding reproducibility task. This inventory does not claim those interactive functions already exist as repository files.

## 12. Version-control checkpoints

Representative completed checkpoints include discovery methods (`d81d221`), discovery results (`70dbeeb`), annotation (`df2052e`), binding audits (`0e632d6`, `1e97c4e`, `e68ff2c`, `2598c9c`), external expression audit (`fa0bd01`), control structure/QC (`de5e958`), first control counts (`dcb9a99`) and second control counts (`87f461d`). The second-control commit was independently checked through GitHub during this session.

When a push was rejected because origin/main had newer commits, I created a backup branch, fetched and rebased, then pushed successfully. Whitespace checks were used before commits. Generated text whitespace was cleaned without changing numeric results.

## 13. Remaining work and interpretation boundaries

1. Confirm SRR29634180 alignment completion; assess mapping, strand assignment and counts.
2. Process the other five knockdowns sequentially and review all sample QC.
3. Save the reusable validation-processing code, commands and software versions.
4. Freeze discovery-only signature selection, direction, weighting/scoring and coverage rules before evaluating external effects. Record exclusions such as whether SREBF2 itself is allowed in the score and define comparator benchmarks.
5. Assemble the eight-sample validation matrix; inspect library sizes, sample relationships and treatment/technical metadata.
6. Apply the frozen score, report each hairpin comparison with uncertainty, and account for shared controls and the small sample size.
7. Assess response concordance and sensitivity without retuning the signature to validation performance.
8. Add binding support only if evidence quality and coordinate compatibility justify it; use wording consistent with the actual evidence.
9. Prepare figures, reproducible methods and a manuscript after assessing transferability and novelty.

No external score-performance result, validated direct-target list or publication-ready final claim has been established yet. Completed processing provides the foundation for those tests. Latest reported resources were 15 GiB total RAM, 9.6 GiB available and 387 GB free disk; these were a snapshot before further extraction, not a current live measurement.

## Groupmate summary

I have built and checked the human reference, processed all 12 discovery RNA-seq samples, validated a 78,932-gene count matrix, completed DESeq2 primary and sensitivity analyses, annotated all tested genes and audited candidate binding and external-expression datasets. Both external HepG2 controls are fully aligned and counted. The first knockdown has passed file and read QC and is ready for alignment assessment. The next deliverables are complete external counts, a frozen discovery-derived signature and an honest evaluation of its external performance.
