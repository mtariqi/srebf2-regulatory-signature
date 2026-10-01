# SREBF2 regulatory signature: rationale, methods, and reproducibility report

**Author/project owner:** Tariqul Islam
**Checkpoint date:** September 30, 2026
**Repository:** https://github.com/mtariqi/srebf2-regulatory-signature
**Scope:** Computational research using real public data; intended duration 8â€“10 weeks.

## 1. Executive assessment

The project has completed repository initialization, a first dataset audit, selection of 12 human discovery runs, environment capture, and extraction and initial quality assessment of one paired-end pilot run, SRR28966297. This report documents the rationale and methods through that checkpoint. It is not a completed biological study.

The pilot generated 28,254,500 spots and 56,509,000 reads, with 28,254,500 reads in each mate file. Both mates have 150-base reads and 50% GC. FastQC and MultiQC finished with exit status 0. A provisional decision was made to retain the reads untrimmed for pilot alignment, while retaining the tile-quality warning and duplication flags for later assessment. Full mate-identifier validation passed for all 28,254,500 read pairs. No reference download, alignment, gene counting, differential expression, binding integration, signature construction, or independent validation has been completed in the visible session.

Evidence in this report comes primarily from terminal outputs and metadata supplied by the project owner. The report author has not independently inspected the full local SRA archive, FASTQ files, environment lockfiles, or complete repository at this checkpoint. Actual hash values and installed package builds should be read from the recorded files, not inferred from this document.

## 2. Scientific rationale

SREBF2 is the gene encoding SREBP2. The project focuses on an evidence-supported transcriptional response associated with this regulator. RNA expression of the regulator and expression of its response genes represent different measurements: the former measures the transcript encoding the regulator, whereas the latter summarizes downstream transcriptional behavior. A target-based score may therefore provide complementary information, but that possibility requires independent validation.

**Research question:** Can a binding-supported, perturbation-derived transcriptional signature capture SREBF2-associated responses more consistently across independent datasets than SREBF2 RNA expression alone?

**Working hypothesis:** A signed score constructed from reproducible SREBF2-responsive genes with appropriate binding support will distinguish relevant perturbation conditions and transfer across independent studies better than, or complement, SREBF2 transcript abundance. This is a hypothesis, not an observed result.

**Candidate research gap:** The cross-study reproducibility and transferability of a frozen SREBF2 signature need evaluation. A systematic recent-literature review remains necessary before claiming this is a novel gap. This report does not establish novelty or provide a comprehensive review of recent SREBF2 research.

**Computational approach:** Use perturbation RNA-seq to nominate responsive genes; integrate existing binding evidence; freeze the gene list, directions, and scoring method using discovery data only; evaluate independent validation data without refitting. Compare the frozen score with SREBF2 expression and a prespecified general cholesterol-related signature. SREBF1 perturbation samples may later assess specificity, but are excluded from the primary 12-run discovery subset.

**Interpretation limits:** Binding plus perturbation support identifies candidate regulatory targets. It does not by itself establish direct regulation, causality, nuclear SREBP2 activation, cholesterol concentration, or clinical utility. No wet-lab validation is included. Synthetic data may be used for software testing only, never to generate the project's biological findings. Wnt data are outside this project.

## 3. Dataset audit and discovery design

### 3.1 Public accession and inputs

- GEO series: **GSE267018**.
- BioProject: **PRJNA1109350**.
- SRA study: **SRP506497**.
- GEO SOFT metadata: `metadata/source/GSE267018_family.soft.gz`.
- Run Selector export: `metadata/GSE267018_SraRunTable.csv`.
- Audit record: `metadata/GSE267018_protocol_audit.txt`.
- Selected run manifest: `metadata/GSE267018_primary_runs.csv`.

The supplied metadata show 30 runs: 18 human HeLa samples and 12 mouse liver samples. Human samples comprise wild type, SREBP2 knockout, and SREBP1 knockout under FBS and cholesterol-depletion conditions, with three labeled replicates per combination. Mouse samples comprise wild type and SREBP2 F167A/F167A liver under fasting and refeeding conditions. Species and experimental contexts must remain separate in primary analyses.

The audit reported poly-A enrichment and paired-end 150-base Illumina sequencing. The deposited study describes HISAT2 and featureCounts processing. Exact annotation release, library strandedness, treatment duration, batch structure, and independence of biological replicates remain unresolved. Replicate labels and distinct BioSamples alone do not prove independent biological replication.

The SOFT record lists HeLa and liver FPKM matrices, with sample-level supplementary files marked NONE. Raw reads are available through SRA. FPKM matrices are unsuitable as direct DESeq2 count input; the project is therefore acquiring raw reads to derive gene-level counts. Availability of another suitable count source has not been exhaustively excluded.

### 3.2 Selected primary runs

| Run | GEO sample | Genotype | Treatment |
|---|---|---|---|
| SRR28966285 | GSM8258553 | SREBP2-/- | Cholesterol depletion |
| SRR28966286 | GSM8258552 | SREBP2-/- | Cholesterol depletion |
| SRR28966287 | GSM8258554 | SREBP2-/- | Cholesterol depletion |
| SRR28966289 | GSM8258551 | SREBP2-/- | FBS |
| SRR28966290 | GSM8258550 | SREBP2-/- | FBS |
| SRR28966291 | GSM8258548 | Wild type | Cholesterol depletion |
| SRR28966292 | GSM8258546 | Wild type | Cholesterol depletion |
| SRR28966293 | GSM8258545 | Wild type | FBS |
| SRR28966294 | GSM8258547 | Wild type | Cholesterol depletion |
| SRR28966295 | GSM8258549 | SREBP2-/- | FBS |
| SRR28966296 | GSM8258544 | Wild type | FBS |
| SRR28966297 | GSM8258543 | Wild type | FBS |

The selected archive sizes sum to approximately 31.78 GB according to the Run Selector Bytes field. That is not a forecast of total working disk usage: uncompressed FASTQ, temporary extraction files, alignments, and references require additional space.

The provisional primary comparison is SREBF2 knockout versus wild type under cholesterol depletion. The FBS comparison and genotype-by-treatment interaction are secondary analyses. A prospective factorial model is `~ genotype + treatment + genotype:treatment`, subject to metadata and design-rank checks. This model has not yet been fitted. Batch terms should be added only when documented and estimable. Replicate numbers must not automatically be treated as paired blocks.

### 3.3 Other candidate datasets

GSE267019 and GSE282800 remain binding-evidence candidates; GSE201466 remains an independent-validation candidate; GSE287355 remains a context-transfer candidate requiring subseries resolution. Their assay suitability, species, sample independence, perturbation identity, and overlap with discovery data have not been established. They are leads, not approved inputs.

## 4. Repository and computational environment

### 4.1 Dedicated repository

The working folder is `/home/mtariq/srebf2-regulatory-signature`, and the user operates fish shell. An initial folder inherited a Git repository rooted at `/home/mtariq`. A dedicated repository was then initialized inside the project folder, and its root was verified. The parent repository was not deleted. Do not stage the user's home directory or unrelated files.

Confirmed user-reported commits:

| Commit | Meaning |
|---|---|
| eda8685 | Initial real-data research plan and tracking |
| d7f1f2b | Metadata-source checkpoint, as recorded in session history |
| 75eac7d | RNA-seq environment capture |
| 70a1d6a | Metadata audit and 12-run human discovery selection |

The pilot report, environment records, checksums, QC summary, and pair validation were committed and pushed in commit 36a8f6b. No direct GitHub write by the assistant is claimed.

Initial tracked material comprises README, research plan, progress tracker, decision log, results index, dataset manifest, sample metadata template, and `.gitignore`. TSV line endings were normalized from CRLF to LF after `git diff --cached --check` flagged carriage returns. New CSV/TSV writers should explicitly use `lineterminator="\n"`.

### 4.2 Hardware and resource choices

The user reported 531 GB available disk space, 15 GiB total RAM, approximately 4.9 GiB available RAM, eight logical processors, and no swap at the initial resource check. Resources fluctuate. Four threads were used for extraction and compression; two FastQC workers processed the mates. Samples will be processed sequentially initially. HISAT2, samtools, and featureCounts are the planned alignment/counting tools; their presence in the environment does not mean those steps have run.

### 4.3 Environments

`srebf2-rnaseq` contains Python 3.12 and FastQC, MultiQC, HISAT2, samtools, subread/featureCounts, SRA tools, and pigz. Versions recorded in the session include FastQC 0.12.1, MultiQC 1.35, HISAT2 2.2.3, samtools 1.24, subread 2.1.1, SRA tools 3.4.1, and pigz 2.8. The committed environment files are authoritative for exact package versions/builds.

A separate environment, `srebf2-sra-test`, was created with SRA tools 3.2.1 during troubleshooting. Successful pilot extraction used that environment and the correct explicit archive path. The original 3.4.1 extraction failure also used an incorrect outer directory, so the observations cannot isolate a version defect. Version 3.2.1 is a working extraction configuration, not proof that 3.4.1 cannot process the archive.

Portable version-oriented exports and platform-specific explicit exports were recorded. The latter are strongest for replay on Linux with the same package platform. Archive URLs or old builds may eventually become unavailable; retaining exports, checksums, logs, and source accessions remains important.

## 5. Executed pilot methods and observations

### 5.1 Download and validation

SRR28966297 is the wild-type/FBS pilot, GEO sample GSM8258543. `prefetch` 3.4.1 reported successful HTTPS download with full base-quality preference and zero unresolved dependencies. The archive is approximately 2.72 GB (displayed as 2.6G by `ls -lh`). `vdb-validate` reported metadata and column MD5 checks passing and a consistent database.

The actual archive location was:

`data/raw/sra/SRR28966297/SRR28966297/SRR28966297.sra`

This location must be discovered rather than assumed. The observed nesting differs from the intended layout.

### 5.2 Troubleshooting audit

1. Extraction with 3.4.1 and the outer directory reported a missing QUALITY column.
2. Diagnostics then used a presumed `.sra` path that did not exist; `vdb-dump` crashed during one diagnostic invocation.
3. An isolated 3.2.1 environment was installed. Inspection of the nonexistent path returned status 0 with no meaningful output; this was not evidence of successful archive inspection.
4. Extraction with 3.2.1 and that nonexistent path failed with resolution/404 errors.
5. A filesystem listing located the nested archive.
6. Extraction using 3.2.1 and the correct absolute file path succeeded.

The correction was both a valid file path and a different environment. The root cause of the initial QUALITY-column message remains unresolved. A reported upstream 3.4.1 crash is contextual evidence only, not a diagnosis of this case.

### 5.3 Extraction, compression, and provenance

`fasterq-dump --split-files --threads 4` produced two approximately 12G uncompressed files. Reported counts were:

| Quantity | Value |
|---|---:|
| Spots read | 28,254,500 |
| Reads read | 56,509,000 |
| Reads written | 56,509,000 |
| Extraction exit status | 0 |

`pigz -p 4` was supplied for compression, followed by `pigz -t` integrity testing. The compressed files subsequently existed and were successfully processed by FastQC. A separate integrity-test exit status was not pasted. SHA-256 values were written to `metadata/SRR28966297_fastq.sha256`; the actual values were not shown in the conversation and are not invented here.

### 5.4 FastQC and MultiQC

FastQC analyzed both compressed mates with two threads. MultiQC 1.35 found two FastQC reports and generated `results/qc/pilot/multiqc/multiqc_report.html`. Both tools returned exit status 0. Tool execution success is distinct from a passing biological/technical QC assessment.

| Measure or module | Mate 1 | Mate 2 |
|---|---:|---:|
| Total sequences | 28,254,500 | 28,254,500 |
| Read length | 150 | 150 |
| GC | 50% | 50% |
| Sequences flagged poor quality | 0 | 0 |
| Per-base quality | PASS | PASS |
| Per-sequence quality | PASS | PASS |
| Per-base N content | PASS | PASS |
| Per-sequence GC | PASS | PASS |
| Overrepresented sequences | PASS | PASS |
| Adapter module | PASS | PASS |
| Per-base composition | FAIL | FAIL |
| Duplication module | FAIL | FAIL |
| Tile quality | WARN | PASS |
| FastQC estimated deduplicated percentage | 41.9317% | 42.6069% |
| Maximum universal adapter signal | 0.0168% | 0.0223% |
| Maximum small-RNA 3â€² adapter signal | 0.0013% | 0.0015% |
| Maximum small-RNA 5â€² adapter signal | 0.0001% | 0.0001% |
| Maximum Nextera signal | 0.0009% | 0.0010% |
| Maximum PolyA signal | 0.5951% | 0.5986% |
| Maximum PolyG signal | 0.0357% | 0.0388% |
| Worst tile deviation | -6.8913, tile 2606, bases 40â€“44 | -2.1036, tile 2656, base 2 |

Base composition is strongly uneven at the first nine positions, approaches balance at bases 10â€“19, and is close to 25% per base by positions 20â€“24. This pattern is consistent with RNA-seq priming bias, although the QC pattern alone does not prove its mechanism. FastQC's duplication estimate is sequence-based and cannot distinguish all biological abundance effects from PCR amplification. It should not be treated as a measured PCR-duplicate fraction.

The tile statistic is a quality deviation relative to the average at the same read position, not an absolute Phred score. The mate-1 warning therefore does not imply Phred 6.89 reads. Neither the number of affected reads nor their absolute quality distribution has been quantified here.

**Provisional processing decision:** retain the untrimmed pilot reads and do not remove duplicates solely to satisfy FastQC. Reassess after alignment, including mapping rate, concordant pairs, assignment rate, strandedness, and unexpected contamination. This is a pilot decision; other samples require their own QC.

## 6. Reproduction commands

The following blocks are fish-compatible. Run from the dedicated project root. Stop when any command fails and inspect its log. Historical failed commands are described above; these blocks consolidate the successful approach and avoid assuming archive nesting. Do not rerun extraction over existing outputs or recreate an existing environment without reviewing its state.

### 6.1 Repository setup and safeguards

For a fresh machine, clone the existing repository rather than initializing over another repository:

```fish
git clone https://github.com/mtariqi/srebf2-regulatory-signature.git ~/srebf2-regulatory-signature
cd ~/srebf2-regulatory-signature
git rev-parse --show-toplevel
git remote -v
git status --short
```

The historical dedicated initialization was `git init -b main`, followed by adding the GitHub origin. Do not repeat initialization or remote creation inside an existing clone.

Relevant ignore rules:

```gitignore
.venv/
__pycache__/
data/raw/
data/processed/
*.fastq
*.fastq.gz
*.fq.gz
*.bam
*.bai
*.cram
work/
.nextflow*
```

### 6.2 Environment creation and capture

```fish
conda create -n srebf2-rnaseq --override-channels -c conda-forge -c bioconda --strict-channel-priority python=3.12 fastqc multiqc hisat2 samtools subread sra-tools pigz
conda activate srebf2-rnaseq
mkdir -p environment
conda env export --no-builds | string match -v 'prefix:*' > environment/rnaseq.yml
conda list --explicit > environment/rnaseq-linux-explicit.txt

conda create -n srebf2-sra-test --override-channels -c conda-forge -c bioconda --strict-channel-priority sra-tools=3.2.1
conda env export -n srebf2-sra-test --no-builds | string match -v 'prefix:*' > environment/sra-extraction.yml
conda list -n srebf2-sra-test --explicit > environment/sra-extraction-linux-explicit.txt
```

The creation commands above describe the historical setup; an unpinned solver may produce a different environment later. For replay, use the captured exports, choosing either the YAML or explicit file per environment:

```fish
conda env create -n srebf2-rnaseq-replay -f environment/rnaseq.yml
conda env create -n srebf2-sra-replay -f environment/sra-extraction.yml
```

Alternatively, on compatible Linux:

```fish
conda create -n srebf2-rnaseq-replay --file environment/rnaseq-linux-explicit.txt
conda create -n srebf2-sra-replay --file environment/sra-extraction-linux-explicit.txt
```

### 6.3 Source metadata and selection

```fish
mkdir -p metadata/source
curl --fail --location --retry 3 --output metadata/source/GSE267018_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267018/soft/GSE267018_family.soft.gz
sha256sum metadata/source/GSE267018_family.soft.gz > metadata/source/GSE267018_family.soft.gz.sha256
gzip -dc metadata/source/GSE267018_family.soft.gz | rg '^\^SAMPLE|^!Sample_title|^!Sample_organism_ch1|^!Sample_characteristics_ch1|^!Sample_supplementary_file|^!Series_supplementary_file'
```

The Run Selector table was obtained manually by searching PRJNA1109350, clearing sample filters, and exporting metadata for all 30 runs. Preserve that export; the study may acquire updated records later. The local copy command was:

```fish
cp ~/Downloads/SraRunTable.csv metadata/GSE267018_SraRunTable.csv
```

This reconstructed equivalent selection code reproduces the documented inclusion rule. It is not asserted to be the byte-for-byte original selection script. Use a separate output first to compare against the committed manifest.

```fish
python3 -c '
import csv
from pathlib import Path
p = Path("metadata/GSE267018_SraRunTable.csv")
with p.open(newline="") as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames
    rows = list(reader)
selected = [r for r in rows if r["Organism"] == "Homo sapiens" and r["genotype"] in {"Wild type", "SREBP2-/-"}]
assert len(selected) == 12, f"Expected 12 runs; found {len(selected)}"
assert len({r["Run"] for r in selected}) == 12, "Duplicate runs"
out = Path("metadata/GSE267018_primary_runs_reproduced.csv")
with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected)
print("Selected runs:", len(selected))
print("Reported archive size (GB):", round(sum(int(r["Bytes"]) for r in selected) / 1e9, 2))
for r in selected:
    print(r["Run"], r["Sample Name"], r["genotype"], r["treatment"])
'
```

### 6.4 Pilot download, path discovery, validation, and extraction

Historical download command:

```fish
conda activate srebf2-rnaseq
mkdir -p data/raw/sra data/raw/fastq data/raw/tmp logs results/qc/pilot
prefetch SRR28966297 -O data/raw/sra/SRR28966297 --max-size 10G > logs/SRR28966297_prefetch.log 2>&1
set download_status $status
echo "Download exit status: $download_status"
tail -n 20 logs/SRR28966297_prefetch.log
```

Discover the actual local file before proceeding. This path-discovery safeguard is a consolidation of the troubleshooting steps:

```fish
set sra_files (find "$PWD/data/raw/sra" -type f -name 'SRR28966297.sra')
count $sra_files
# Continue only if exactly one archive was found.
set sra_file $sra_files[1]
ls -lh "$sra_file"
vdb-validate "$sra_file" > logs/SRR28966297_validate.log 2>&1
set validate_status $status
echo "Validation exit status: $validate_status"
tail -n 20 logs/SRR28966297_validate.log
```

After validation succeeds:

```fish
conda activate srebf2-sra-test
fasterq-dump --version
fasterq-dump "$sra_file" --split-files --threads 4 --outdir data/raw/fastq --temp data/raw/tmp > logs/SRR28966297_extract_sra321.log 2>&1
set extract_status $status
echo "Extraction exit status: $extract_status"
tail -n 20 logs/SRR28966297_extract_sra321.log
```

The variables persist within the same fish session. In a new shell, set `sra_file` again. The successful historical path was `$PWD/data/raw/sra/SRR28966297/SRR28966297/SRR28966297.sra`.

### 6.5 Compression, integrity, and checksums

```fish
conda activate srebf2-rnaseq
pigz -p 4 data/raw/fastq/SRR28966297_1.fastq data/raw/fastq/SRR28966297_2.fastq
and pigz -t data/raw/fastq/SRR28966297_1.fastq.gz data/raw/fastq/SRR28966297_2.fastq.gz
sha256sum data/raw/fastq/SRR28966297_1.fastq.gz data/raw/fastq/SRR28966297_2.fastq.gz > metadata/SRR28966297_fastq.sha256
```

Compression replaces the original FASTQ files with gzip files. A hash permits future identity checks but does not itself verify biological provenance. Independently recompressed files may have different gzip bytes despite identical decompressed reads; compare decompressed content if necessary.

### 6.6 FastQC and MultiQC

```fish
fastqc --threads 2 --outdir results/qc/pilot data/raw/fastq/SRR28966297_1.fastq.gz data/raw/fastq/SRR28966297_2.fastq.gz > logs/SRR28966297_fastqc.log 2>&1
set qc_status $status
echo "FastQC exit status: $qc_status"
tail -n 20 logs/SRR28966297_fastqc.log
# Continue only when FastQC succeeds.
multiqc results/qc/pilot --outdir results/qc/pilot/multiqc > logs/SRR28966297_multiqc.log 2>&1
set multiqc_status $status
echo "MultiQC exit status: $multiqc_status"
tail -n 20 logs/SRR28966297_multiqc.log
```

## 7. QC parsing code used in the session

Save the following blocks as Python scripts under `scripts/` and run them from the project root. These consolidate the executed inline Python commands with small readability changes. No third-party Python dependencies are required.

### 7.1 Basic summary: `scripts/pilot_qc_summary.py`

```python
from pathlib import Path
from zipfile import ZipFile

reports = sorted(Path("results/qc/pilot").glob("*_fastqc.zip"))
assert len(reports) == 2, f"Expected 2 FastQC archives; found {len(reports)}"
output = []
for path in reports:
    with ZipFile(path) as z:
        summary = next(n for n in z.namelist() if n.endswith("/summary.txt"))
        data = next(n for n in z.namelist() if n.endswith("/fastqc_data.txt"))
        output.append(f"\n## {path.name}\n")
        output.append(z.read(summary).decode())
        output.append("\nBasic statistics:\n")
        for line in z.read(data).decode().splitlines():
            if line.startswith(("Total Sequences\t", "Sequence length\t", "%GC\t", "Sequences flagged as poor quality\t")):
                output.append(line + "\n")
text = "".join(output)
Path("results/qc/pilot/QC_SUMMARY.txt").write_text(text)
print(text)
```

### 7.2 Full flag extraction: `scripts/pilot_qc_flag_details.py`

```python
from pathlib import Path
from zipfile import ZipFile

wanted = {"Per base sequence content", "Sequence Duplication Levels", "Per tile sequence quality", "Adapter Content"}
output = []
for path in sorted(Path("results/qc/pilot").glob("*_fastqc.zip")):
    output.append(f"\n## {path.name}\n")
    with ZipFile(path) as z:
        name = next(n for n in z.namelist() if n.endswith("/fastqc_data.txt"))
        keep = False
        for line in z.read(name).decode().splitlines():
            if line.startswith(">>") and line != ">>END_MODULE":
                keep = line[2:].split("\t")[0] in wanted
            if keep:
                output.append(line + "\n")
            if line == ">>END_MODULE":
                keep = False
text = "".join(output)
Path("results/qc/pilot/QC_FLAG_DETAILS.txt").write_text(text)
print(text)
```

### 7.3 Compact flag metrics: `scripts/pilot_qc_compact.py`

```python
from pathlib import Path
from zipfile import ZipFile

for p in sorted(Path("results/qc/pilot").glob("*_fastqc.zip")):
    print("\n" + p.name)
    with ZipFile(p) as z:
        name = next(n for n in z.namelist() if n.endswith("/fastqc_data.txt"))
        text = z.read(name).decode()
    for block in text.split(">>END_MODULE"):
        lines = block.strip().splitlines()
        if not lines:
            continue
        title = lines[0]
        rows = [x.split("\t") for x in lines[1:] if x and not x.startswith("#")]
        if title.startswith(">>Per tile sequence quality") and rows:
            worst = min(rows, key=lambda r: float(r[2]))
            print(title)
            print("Worst tile/base/quality deviation:", *worst)
        elif title.startswith(">>Sequence Duplication Levels"):
            print(title)
            for line in lines:
                if line.startswith("#Total"):
                    print(line)
        elif title.startswith(">>Adapter Content") and rows:
            headers = next(x for x in lines if x.startswith("#Position")).split("\t")
            for i, label in enumerate(headers[1:], 1):
                print(f"Maximum {label}: {max(float(r[i]) for r in rows):.4f}%")
        elif title.startswith(">>Per base sequence content"):
            print(title)
            print("Base / G / A / T / C:")
            for row in rows[:12]:
                print(*row)
```

### 7.4 Mate validation: `scripts/validate_pilot_pairs.py` â€” completed successfully

```python
import gzip
from pathlib import Path

root = Path("data/raw/fastq")
count = 0

def read_record(handle):
    header = handle.readline()
    if not header:
        return None
    seq, plus, qual = (handle.readline() for _ in range(3))
    if not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError("Invalid FASTQ record")
    if not seq or not qual or len(seq.rstrip()) != len(qual.rstrip()):
        raise ValueError("Incomplete record or sequence/quality length mismatch")
    identifier = header.split()[0]
    if identifier.endswith((b"/1", b"/2")):
        identifier = identifier[:-2]
    return identifier

with gzip.open(root / "SRR28966297_1.fastq.gz", "rb") as a, gzip.open(root / "SRR28966297_2.fastq.gz", "rb") as b:
    while True:
        r1, r2 = read_record(a), read_record(b)
        if r1 is None and r2 is None:
            break
        if r1 is None or r2 is None or r1 != r2:
            raise ValueError(f"Mate mismatch at pair {count + 1}: {r1!r}, {r2!r}")
        count += 1
assert count == 28254500, f"Unexpected pair count: {count}"
message = f"PASS: {count:,} matched read pairs; FASTQ records valid.\n"
Path("results/qc/pilot/PAIR_VALIDATION.txt").write_text(message)
print(message)
```

This streams the files with constant memory and checks structural FASTQ validity, sequence/quality length equality, matching first-token identifiers, and expected count. It is not a comprehensive character/alphabet validator. Matching identifiers at every position were confirmed for all 28,254,500 pairs.

## 8. Outputs and traceability

| Artifact | Purpose | Status at this checkpoint |
|---|---|---|
| SOFT archive and `.sha256` | Source metadata and identity | Recorded previously |
| SRA Run Selector CSV | Original run metadata | Committed |
| Primary run CSV | Inclusion record | Committed |
| Protocol audit | Dataset questions and protocol evidence | Committed |
| RNA-seq YAML and explicit export | Processing environment | Committed |
| Extraction YAML and explicit export | Working extraction environment | Committed in 36a8f6b |
| Pilot SRA archive | Downloaded source reads | Validated locally |
| Pilot compressed FASTQ mates | Extracted reads | FastQC processed |
| FASTQ `.sha256` file | Input identity | Written; values not pasted |
| Prefetch/validation/extraction logs | Acquisition evidence | User-reported outputs |
| Two FastQC ZIP/HTML reports | Raw-read QC | Generated |
| MultiQC report/data | Aggregated QC | Generated |
| QC_SUMMARY.txt | Human-readable summary | Generated |
| QC_FLAG_DETAILS.txt | Detailed module data | Generated |
| PAIR_VALIDATION.txt | Full mate-identifier check | PASS; committed in 36a8f6b |

Raw reads and large intermediate files must remain outside Git. Retain accessions, checksums, scripts, environment exports, logs, and compact QC evidence in version control. For every later result, record the input hashes, reference/annotation release, command, software version, generating Git commit, QC status, and interpretation. Avoid claiming provenance solely from a filename.

Paired-read integrity validation. Both compressed FASTQ files for SRR28966297 were streamed using Python’s gzip module. The validation checked FASTQ header and separator structure, sequence–quality length equality, mate-identifier agreement, and the expected number of pairs. All 28,254,500 read pairs passed, with no detected mate mismatches or incomplete records.
## 9. Pending work and prespecified analysis safeguards

1. Completed: pair validation passed and the pilot checkpoint was pushed.
2. Download a matched GRCh38 genome and annotation; freeze exact release, chromosome naming, source URLs, and hashes.
3. Establish HISAT2 index strategy appropriate to available RAM and document splice-site use.
4. Align the pilot; inspect mapping, concordant pairing, multimapping, gene assignment, and strand evidence before selecting featureCounts strandedness.
5. Acquire and process the remaining 11 primary runs using the same frozen workflow; review sample-specific QC and contamination.
6. Confirm biological replication, treatments, batch variables, and count-column/sample mapping before DESeq2.
7. Perform discovery differential expression on raw integer gene counts. Report effect sizes, uncertainty, and multiple-testing-adjusted significance. Freeze contrast definitions and filtering rules before selecting signature genes.
8. Audit binding datasets, genome coordinates, controls, and study overlap. Binding integration has not yet been implemented.
9. Exclude SREBF2 itself from the target-based score; freeze discovery gene membership, direction, scoring and missing-gene handling before validation.
10. Evaluate independent validation without refitting; compare with regulator RNA and a prespecified cholesterol comparator; assess specificity and robustness.
11. Generate report figures and a reproducible release. Negative or inconsistent findings remain reportable results.

Eight-week sequence: literature/data audit; acquisition/QC; discovery counts/DE; binding/signature; independent validation; comparator analyses; robustness; report/release. Weeks 9â€“10 provide processing and validation buffer. No completed result should be inferred from this schedule.

## 10. Methods text suitable for a later manuscript

> At the initial processing checkpoint, public metadata for GEO series GSE267018 and BioProject PRJNA1109350 were audited. Twelve human HeLa RNA-seq runs representing wild-type and SREBP2-knockout cells under FBS and cholesterol-depletion conditions were selected, with three labeled replicates per condition. Biological replication and additional design variables remained under audit. Raw reads for the wild-type/FBS pilot SRR28966297 were downloaded using SRA Toolkit 3.4.1 and passed archive validation. FASTQ extraction succeeded using SRA Toolkit 3.2.1 with the correct explicit archive path, yielding 28,254,500 paired spots. Reads were compressed using pigz, and SHA-256 checksums were recorded. Raw-read quality was assessed using FastQC and summarized using MultiQC 1.35. Both mates passed sequence-quality and adapter modules, while base-composition and duplication modules failed; mate 1 also showed a tile-quality warning. The reads were provisionally retained untrimmed for pilot alignment. No alignment or differential-expression results were available at this checkpoint.

Before manuscript submission, replace unresolved version/design details with verified evidence and update this paragraph with completed methods. Do not present this preliminary paragraph as the methods for a finished signature study.

## 11. Primary-source references (APA style)

Andrews, S. (n.d.). *FastQC: A quality control tool for high throughput sequence data*. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Babraham Bioinformatics. (n.d.). *Duplicate sequences*. FastQC documentation. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html

Babraham Bioinformatics. (n.d.). *Per base sequence content*. FastQC documentation. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html

National Center for Biotechnology Information. (n.d.). *GSE267018* [Data set metadata]. Gene Expression Omnibus. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267018

National Center for Biotechnology Information. (n.d.). *08. Prefetch and fasterq dump*. SRA Tools documentation. https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump

National Center for Biotechnology Information. (n.d.). *SRA Tools* [Computer software]. GitHub. https://github.com/ncbi/sra-tools

Dates and original study authors were not inferred where unverified. The GEO reference identifies the metadata record; the original study paper and complete dataset citation must be added after bibliographic verification. Documentation was consulted on September 30, 2026. This reference list supports the completed processing workflow; it is not a recent-literature review establishing the project's novelty.
