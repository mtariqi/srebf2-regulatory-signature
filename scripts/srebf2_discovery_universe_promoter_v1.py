#!/usr/bin/env python3
"""Extract SREBP2 promoter bigWig signal for the frozen GSE267018 evaluated universe.
Run from the root of srebf2-regulatory-signature with srebf2-binding environment.
"""
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path

import pyBigWig

AUDIT = Path('results/signature/discovery_v1/SIGNATURE_SELECTION_AUDIT.tsv')
SIGNATURE = Path('results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv')
GTF = Path('data/raw/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf.gz')
BW_DIR = Path('data/processed/binding/GSE282800')
OUT_DIR = Path('results/binding/GSE282800')
OUT = OUT_DIR / 'SREBF2_DISCOVERY_UNIVERSE_PROMOTER_SIGNAL_v1.tsv'
MANIFEST = OUT_DIR / 'DISCOVERY_UNIVERSE_PROMOTER_MANIFEST_v1.json'
UPSTREAM, DOWNSTREAM = 2000, 500
PRIMARY = {f'chr{x}' for x in range(1, 23)} | {'chrX', 'chrY'}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load_table(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows:
        raise ValueError(f'Empty table: {path}')
    return rows


def parse_gtf(path, wanted):
    genes = {}
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) != 9 or p[2] != 'gene':
                continue
            match = re.search(r'(?:^|;)\s*gene_id\s+"([^"]+)"', p[8])
            if not match or match.group(1) not in wanted:
                continue
            gene_id = match.group(1)
            if gene_id in genes:
                raise ValueError(f'Duplicate gene-level GTF entry: {gene_id}')
            strand = p[6]
            if strand not in ('+', '-'):
                continue
            tss = int(p[3]) - 1 if strand == '+' else int(p[4]) - 1
            genes[gene_id] = (p[0], strand, tss)
    return genes


def score(bw, chrom, start, end):
    return (bw.stats(chrom, start, end, type='mean', exact=True)[0],
            bw.stats(chrom, start, end, type='max', exact=True)[0])


def as_text(value):
    return '' if value is None else str(value)


def main():
    audit = load_table(AUDIT)
    signature = load_table(SIGNATURE)
    required = {'gene_id', 'gene_symbol', 'gene_type', 'baseMean',
                'signature_selected', 'direction', 'log2FC_primary'}
    missing = required - audit[0].keys()
    if missing:
        raise ValueError(f'Missing audit columns: {sorted(missing)}')
    ids = [r['gene_id'] for r in audit]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate gene_id values in audit')
    if len(audit) != 17158 or len(signature) != 1326:
        raise ValueError(f'Unexpected frozen counts: audit={len(audit)}, signature={len(signature)}')
    sig_ids = {r['gene_id'] for r in signature}
    selected_ids = {r['gene_id'] for r in audit if r['signature_selected'].strip().lower() == 'true'}
    if sig_ids != selected_ids:
        raise ValueError('Selected audit genes do not match frozen signature IDs')
    gtf = GTF if GTF.exists() else GTF.with_suffix('')
    if not gtf.exists():
        raise FileNotFoundError('GENCODE v47 GTF not found')
    genes = parse_gtf(gtf, set(ids))
    fbs_files = list(BW_DIR.glob('*FBS-IP*.bw'))
    dep_files = list(BW_DIR.glob('*LPPS-IP*.bw'))
    if len(fbs_files) != 1 or len(dep_files) != 1:
        raise ValueError('Expected exactly one FBS-IP and one LPPS-IP bigWig')
    fbs_path, dep_path = fbs_files[0], dep_files[0]
    extra_fields = ['chrom', 'strand', 'tss_0based', 'promoter_start_0based',
                    'promoter_end_0based_exclusive', 'FBS_mean', 'FBS_max',
                    'depletion_mean', 'depletion_max', 'depletion_minus_FBS',
                    'mapping_status', 'primary_nuclear_chromosome',
                    'nonzero_in_either_track']
    out_fields = list(audit[0].keys()) + [f for f in extra_fields if f not in audit[0]]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    statuses = Counter()
    with pyBigWig.open(str(fbs_path)) as fbs, pyBigWig.open(str(dep_path)) as dep:
        if fbs.chroms() != dep.chroms():
            raise ValueError('bigWig chromosome dictionaries differ')
        lengths = fbs.chroms()
        with OUT.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, delimiter='\t', fieldnames=out_fields)
            writer.writeheader()
            for row in audit:
                record = dict(row)
                record.update({f: '' for f in extra_fields})
                record['primary_nuclear_chromosome'] = 'False'
                record['nonzero_in_either_track'] = ''
                annotation = genes.get(row['gene_id'])
                if annotation is None:
                    status = 'unmatched_gene_id'
                else:
                    chrom, strand, tss = annotation
                    record.update(chrom=chrom, strand=strand, tss_0based=tss,
                                  primary_nuclear_chromosome=str(chrom in PRIMARY))
                    if chrom not in lengths:
                        status = 'chromosome_not_in_bigwig'
                    else:
                        start = tss - UPSTREAM if strand == '+' else tss - DOWNSTREAM
                        end = tss + DOWNSTREAM + 1 if strand == '+' else tss + UPSTREAM + 1
                        start, end = max(0, start), min(lengths[chrom], end)
                        record['promoter_start_0based'] = start
                        record['promoter_end_0based_exclusive'] = end
                        if end <= start:
                            status = 'invalid_interval'
                        else:
                            fm, fx = score(fbs, chrom, start, end)
                            dm, dx = score(dep, chrom, start, end)
                            record.update(FBS_mean=as_text(fm), FBS_max=as_text(fx),
                                          depletion_mean=as_text(dm), depletion_max=as_text(dx))
                            if fm is None or dm is None:
                                status = 'mapped_missing_signal'
                            else:
                                record['depletion_minus_FBS'] = dm - fm
                                record['nonzero_in_either_track'] = str(fm > 0 or dm > 0)
                                status = 'mapped_with_signal'
                record['mapping_status'] = status
                statuses[status] += 1
                writer.writerow(record)
    manifest = {
        'analysis': 'discovery_universe_promoter_v1',
        'classification': 'descriptive processed ChIP-seq signal; not peak-level binding',
        'genome': 'GRCh38/hg38', 'annotation': 'GENCODE v47 gene-level TSS',
        'promoter': {'upstream_bp': UPSTREAM, 'downstream_bp': DOWNSTREAM,
                     'coordinate_system': '0-based half-open'},
        'input_rows': len(audit), 'frozen_selected_rows': len(sig_ids),
        'mapping_status_counts': dict(statuses),
        'primary_nuclear_definition': 'chr1-chr22,chrX,chrY',
        'files_sha256': {str(p): sha256(p) for p in (AUDIT, SIGNATURE, gtf, fbs_path, dep_path, OUT)},
        'caveats': ['No replicate-level differential binding inference',
                    'No claim of direct SREBP2 target assignment',
                    'Zero signal is not necessarily evidence of biological absence',
                    'Do not use this result to alter frozen signature membership'],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + '\n')
    print('Evaluated universe:', len(audit))
    print('Frozen signature:', len(sig_ids))
    print('GTF matches:', len(genes))
    print('Mapping statuses:', dict(statuses))
    print('Output:', OUT)
    print('SHA256:', sha256(OUT))
    print('Manifest:', MANIFEST)


if __name__ == '__main__':
    main()
