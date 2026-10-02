"""Run from repository root. Standard-library-only, exact GENCODE ID mapping."""
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

root = Path('results/de/discovery_v1')
gtf = Path('data/raw/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf.gz')
output = root / 'annotated'
if output.exists():
    raise SystemExit('Output exists; preserve it and choose a new output version.')
expected = 'f02ee3e1c8e7fd9be264be6d0b974feb225a1e9d6c81915ee271804b060bd8c0'
def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()
if sha256(gtf) != expected:
    raise SystemExit('GENCODE annotation checksum differs from the recorded reference.')
annotation = {}
with gzip.open(gtf, 'rt') as handle:
    for line in handle:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            raise ValueError('Malformed GTF line')
        if fields[2] != 'gene':
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
        gene = attrs['gene_id']
        if gene in annotation:
            raise ValueError(f'Duplicate GTF gene ID: {gene}')
        annotation[gene] = {'gene_symbol': attrs.get('gene_name', ''),
                            'gene_type': attrs.get('gene_type', ''),
                            'chromosome': fields[0], 'start_1based': fields[3],
                            'end_1based': fields[4], 'strand': fields[6]}
sources = [root / fit / (contrast + '.tsv')
           for fit in ('all12', 'without_SRR28966285')
           for contrast in ('KO_vs_WT_depletion', 'KO_vs_WT_FBS', 'interaction')]
sources += [root / (contrast + '_SENSITIVITY.tsv')
            for contrast in ('KO_vs_WT_depletion', 'KO_vs_WT_FBS', 'interaction')]
for path in sources:
    if not path.is_file():
        raise FileNotFoundError(path)
output.mkdir()
extras = ['gene_symbol', 'gene_type', 'chromosome', 'start_1based', 'end_1based', 'strand', 'annotation_status']
qc = []
srebf2 = []
hashes = {str(gtf): expected}
for source in sources:
    target = output / source.relative_to(root)
    target.parent.mkdir(parents=True, exist_ok=True)
    hashes[str(source)] = sha256(source)
    seen = set()
    mapped = 0
    with source.open(newline='') as incoming, target.open('w', newline='') as outgoing:
        reader = csv.DictReader(incoming, delimiter='\t')
        original = reader.fieldnames
        if not original or 'gene_id' not in original or set(original) & set(extras):
            raise ValueError(f'Unexpected columns: {source}')
        writer = csv.DictWriter(outgoing, fieldnames=original + extras, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in reader:
            gene = row['gene_id']
            if not gene or gene in seen:
                raise ValueError(f'Missing/duplicate gene ID: {source}, {gene}')
            seen.add(gene)
            info = annotation.get(gene)
            mapped += info is not None
            row.update(info or {key: '' for key in extras[:-1]})
            row['annotation_status'] = 'exact_match' if info else 'unmapped'
            writer.writerow(row)
            if row['gene_symbol'] == 'SREBF2':
                srebf2.append({'source': str(source), **row})
    qc.append({'source': str(source), 'rows': len(seen), 'mapped': mapped, 'unmapped': len(seen)-mapped})
    print(f'PASS: {source}; {len(seen):,} rows; {mapped:,} mapped; {len(seen)-mapped} unmapped')
with (output / 'ANNOTATION_QC.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(qc[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(qc)
with (output / 'INPUT_SHA256.json').open('w') as handle:
    json.dump(hashes, handle, indent=2); handle.write('\n')
with (output / 'SREBF2_RESULTS.json').open('w') as handle:
    json.dump(srebf2, handle, indent=2); handle.write('\n')
(output / 'METHODS.txt').write_text(
    'GENCODE v47 primary assembly GTF; compressed reference checksum verified.\n'
    'Only gene features used; exact versioned gene_id matching; no version stripping.\n'
    'Source statistical values and row order preserved; unmapped rows retained and flagged.\n'
    'Duplicate gene IDs rejected; duplicate symbols not collapsed.\n'
    'SREBF2 extraction is descriptive QC, not a target signature or protein-activity measurement.\n')
print('COMPLETE: review ANNOTATION_QC.tsv and SREBF2_RESULTS.json before interpretation.')
