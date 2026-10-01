#!/usr/bin/env python3
"""Run from the project root in the srebf2-rnaseq environment.
Real-data pipeline; one sample at a time. Incomplete work is not silently reused.
"""
import csv
import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path
from zipfile import ZipFile

ROOT = Path.cwd()
REF = Path('data/raw/reference/gencode_v47')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def command(args, log):
    with log.open('w') as f:
        f.write('COMMAND: ' + json.dumps([str(x) for x in args]) + '\n')
        f.flush()
        subprocess.run([str(x) for x in args], stdout=f, stderr=subprocess.STDOUT, check=True)

def pairs(a, b):
    def record(f):
        h = f.readline()
        if not h:
            return None
        s, p, q = f.readline(), f.readline(), f.readline()
        if not h.startswith(b'@') or not p.startswith(b'+') or not s or not q:
            raise ValueError('Invalid FASTQ record')
        if len(s.rstrip(b'\r\n')) != len(q.rstrip(b'\r\n')):
            raise ValueError('Sequence/quality length mismatch')
        name = h.split()[0]
        return name[:-2] if name.endswith((b'/1', b'/2')) else name
    n = 0
    with gzip.open(a, 'rb') as x, gzip.open(b, 'rb') as y:
        while True:
            u, v = record(x), record(y)
            if u is None and v is None:
                break
            if u is None or v is None or u != v:
                raise ValueError(f'Mate mismatch at pair {n + 1}')
            n += 1
    if not n:
        raise ValueError('Empty FASTQs')
    return n

def process(run, gsm):
    qc = Path('results/qc') / run
    complete = qc / 'PIPELINE_COMPLETE.json'
    if complete.exists():
        saved = json.loads(complete.read_text())
        for name, digest in saved['checksums'].items():
            if sha(Path(name)) != digest:
                raise ValueError(f'Completed output changed: {name}')
        subprocess.run(['samtools', 'quickcheck', saved['bam']], check=True)
        print(f'{run}: completed outputs verified; skipping', flush=True)
        return
    reads = [Path(f'data/raw/fastq/{run}_{mate}.fastq.gz') for mate in (1, 2)]
    bam = Path(f'data/processed/alignment/{run}.sorted.bam')
    if qc.exists() or bam.exists() or any(p.exists() or p.with_suffix('').exists() for p in reads):
        raise RuntimeError(f'{run}: incomplete outputs exist. Stop and inspect before resuming; no files were deleted.')
    qc.mkdir(parents=True)
    logs = Path('logs') / run
    logs.mkdir(parents=True, exist_ok=True)
    tmp = Path('data/raw/tmp') / run
    tmp.mkdir(parents=True, exist_ok=True)
    print(f'{run}: downloading', flush=True)
    command(['prefetch', run, '-O', 'data/raw/sra', '--max-size', '10G'], logs / 'prefetch.log')
    archive = (Path('data/raw/sra') / run / f'{run}.sra').resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)
    command(['conda', 'run', '-n', 'srebf2-sra-test', 'vdb-validate', archive], logs / 'validate.log')
    command(['conda', 'run', '-n', 'srebf2-sra-test', 'fasterq-dump', archive,
             '--split-files', '--threads', '4', '--outdir', 'data/raw/fastq', '--temp', tmp], logs / 'extract.log')
    extraction = (logs / 'extract.log').read_text()
    match = re.search(r'spots read\s*:\s*([\d,]+)', extraction)
    if not match:
        raise ValueError('Cannot obtain expected pair count from extraction log')
    expected = int(match.group(1).replace(',', ''))
    command(['pigz', '-p', '4', *[p.with_suffix('') for p in reads]], logs / 'compression.log')
    n = pairs(*reads)
    if n != expected:
        raise ValueError(f'Pair count {n} differs from extraction count {expected}')
    (qc / 'PAIR_VALIDATION.txt').write_text(f'PASS: {n:,} matched read pairs; FASTQ records valid.\n')
    hashes = {str(p): sha(p) for p in reads}
    Path(f'metadata/{run}_fastq.sha256').write_text(''.join(f'{h}  {p}\n' for p, h in hashes.items()))
    command(['fastqc', '--threads', '2', '--outdir', qc, *reads], logs / 'fastqc.log')
    sections = []
    for mate in (1, 2):
        z = qc / f'{run}_{mate}_fastqc.zip'
        with ZipFile(z) as f:
            member = next(x for x in f.namelist() if x.endswith('/summary.txt'))
            sections.append(z.name + '\n' + f.read(member).decode().rstrip())
    (qc / 'QC_SUMMARY.txt').write_text('\n\n'.join(sections) + '\n')
    print(f'{run}: aligning (QC flags remain subject to review)', flush=True)
    align = ['hisat2', '-p', '4', '-x', str(REF / 'hisat2_index/GRCh38_gencode_v47'),
             '--known-splicesite-infile', str(REF / 'gencode.v47.splicesites.tsv'),
             '--rg-id', run, '--rg', f'SM:{gsm}', '--rg', 'PL:ILLUMINA', '-1', str(reads[0]), '-2', str(reads[1])]
    sort = ['samtools', 'sort', '-@', '1', '-m', '768M', '-T', str(tmp / 'sort'), '-o', str(bam), '-']
    with (logs / 'hisat2.log').open('w') as al, (logs / 'sort.log').open('w') as sl:
        p = subprocess.Popen(align, stdout=subprocess.PIPE, stderr=al)
        try:
            q = subprocess.Popen(sort, stdin=p.stdout, stderr=sl)
        except Exception:
            p.stdout.close(); p.terminate(); p.wait(); raise
        p.stdout.close()
        sr, ar = q.wait(), p.wait()
    if ar or sr:
        raise RuntimeError(f'Alignment/sort failed: {ar}/{sr}')
    command(['samtools', 'quickcheck', '-v', bam], logs / 'quickcheck.log')
    command(['samtools', 'index', bam], logs / 'index.log')
    with (qc / f'{run}_flagstat.txt').open('w') as f:
        subprocess.run(['samtools', 'flagstat', '-@', '2', str(bam)], stdout=f, check=True)
    (qc / f'{run}_hisat2_summary.txt').write_text((logs / 'hisat2.log').read_text())
    countdir = qc / 'strand_check'; countdir.mkdir()
    comparison = ['strand_mode\tassigned\ttotal_summary_entries\tassigned_percent']
    assigned0 = None
    for mode in (0, 1, 2):
        out = countdir / f'{run}_s{mode}.txt'
        command(['featureCounts', '-T', '4', '-p', '--countReadPairs', '-B', '-C', '-s', str(mode),
                 '-t', 'exon', '-g', 'gene_id', '-a', REF / 'gencode.v47.primary_assembly.annotation.gtf',
                 '-o', out, bam], logs / f'featureCounts_s{mode}.log')
        with Path(str(out) + '.summary').open() as f:
            rr = list(csv.reader(f, delimiter='\t'))
        vals = {r[0]: int(r[1]) for r in rr[1:]}
        total = sum(vals.values()); assigned = vals['Assigned']
        if not total: raise ValueError('Empty count summary')
        comparison.append(f'{mode}\t{assigned}\t{total}\t{100*assigned/total:.2f}')
        if mode == 0: assigned0 = assigned
    (countdir / 'STRAND_ASSIGNMENT_COMPARISON.tsv').write_text('\n'.join(comparison) + '\n')
    with (countdir / f'{run}_s0.txt').open() as f:
        reader = csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t')
        column = reader.fieldnames[-1]
        counts = [(r['Geneid'], int(r[column])) for r in reader]
    if len(counts) != len(set(g for g, _ in counts)) or any(v < 0 for _, v in counts):
        raise ValueError('Invalid gene IDs or counts')
    if sum(v for _, v in counts) != assigned0: raise ValueError('Assigned count sum mismatch')
    compact = Path(f'results/counts/{run}_gene_counts.tsv')
    with compact.open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_id', run]); writer.writerows(counts)
    (countdir / 'COUNT_VALIDATION.txt').write_text(f'PASS: unique gene IDs and nonnegative integer counts\nGene rows: {len(counts):,}\nGenes with nonzero counts: {sum(v>0 for _,v in counts):,}\nTotal assigned fragments: {assigned0:,}\n')
    hashes[str(bam)] = sha(bam)
    Path(f'metadata/{run}_sorted_bam.sha256').write_text(f'{hashes[str(bam)]}  {bam}\n')
    hashes[str(compact)] = sha(compact)
    complete.write_text(json.dumps({'run': run, 'sample': gsm, 'pairs': n, 'bam': str(bam),
                                   'qc_review': 'pending', 'checksums': hashes,
                                   'alignment_command': align, 'sort_command': sort}, indent=2) + '\n')
    print(f'{run}: processing complete; QC review pending', flush=True)

if __name__ == '__main__':
    for folder in ['data/raw/sra', 'data/raw/fastq', 'data/processed/alignment', 'results/counts', 'metadata', 'logs']:
        Path(folder).mkdir(parents=True, exist_ok=True)
    subprocess.run(['sha256sum', '-c', 'metadata/GRCh38_gencode_v47_hisat2_index.sha256'], check=True)
    subprocess.run(['sha256sum', '-c', 'metadata/gencode_v47_splicesites.sha256'], check=True)
    with Path('metadata/GSE267018_primary_runs.csv').open() as f:
        manifest = {r['Run']: r for r in csv.DictReader(f)}
    queue = Path('metadata/GSE267018_processing_queue.txt').read_text().splitlines()
    if len(queue) != len(set(queue)) or not set(queue).issubset(manifest):
        raise ValueError('Invalid processing queue')
    for run in queue:
        row = manifest[run]
        gsm = row.get('Sample Name') or row.get('Library Name')
        if not gsm or not re.fullmatch(r'GSM\d+', gsm): raise ValueError(f'Invalid GSM for {run}')
        process(run, gsm)
    print('QUEUE COMPLETE. Review all QC before downstream analysis.', flush=True)
