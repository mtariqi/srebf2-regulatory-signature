#!/usr/bin/env python3
"""Freeze a discovery-only, perturbation-derived SREBF2 candidate signature.
Run from the repository root: python scripts/srebf2_build_signature_v1.py
"""
import hashlib
import json
from datetime import date
from pathlib import Path
import pandas as pd

ROOT = Path.cwd()
BASE = ROOT / 'results/de/discovery_v1/annotated'
OUT = ROOT / 'results/signature/discovery_v1'
SOURCES = {
    'primary': BASE / 'all12/KO_vs_WT_depletion.tsv',
    'fbs': BASE / 'all12/KO_vs_WT_FBS.tsv',
    'sensitivity': BASE / 'without_SRR28966285/KO_vs_WT_depletion.tsv',
}
THRESHOLDS = {'primary_fdr_max_exclusive': 0.05, 'primary_abs_log2fc_min_inclusive': 1.0}

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def load(path, tag):
    if not path.is_file():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep='\t', dtype={'gene_id': str, 'gene_symbol': str}, low_memory=False)
    required = {'gene_id', 'log2FoldChange', 'padj', 'gene_symbol', 'gene_type', 'annotation_status'}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}')
    if df.gene_id.isna().any() or df.gene_id.duplicated().any():
        raise ValueError(f'{path}: missing or duplicate gene IDs')
    for c in ['log2FoldChange', 'padj']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df.set_index('gene_id')

if OUT.exists():
    raise SystemExit(f'Refusing to overwrite frozen signature directory: {OUT}')

a = load(SOURCES['primary'], 'primary')
b = load(SOURCES['fbs'], 'fbs')
c = load(SOURCES['sensitivity'], 'sensitivity')
if len(a) != 17158 or set(a.index) != set(b.index) or set(a.index) != set(c.index):
    raise ValueError('Discovery tables must share exactly 17,158 gene IDs')

m = a[['gene_symbol', 'gene_type', 'annotation_status', 'baseMean', 'log2FoldChange', 'padj']].copy()
m = m.rename(columns={'log2FoldChange': 'log2FC_primary', 'padj': 'padj_primary'})
m['log2FC_fbs'] = b.reindex(m.index)['log2FoldChange']
m['padj_fbs'] = b.reindex(m.index)['padj']
m['log2FC_sensitivity'] = c.reindex(m.index)['log2FoldChange']
m['padj_sensitivity'] = c.reindex(m.index)['padj']
m['primary_fdr_pass'] = m.padj_primary.notna() & (m.padj_primary < 0.05)
m['primary_effect_pass'] = m.log2FC_primary.notna() & (m.log2FC_primary.abs() >= 1.0)
m['fbs_direction_agrees'] = (m.log2FC_primary * m.log2FC_fbs > 0).fillna(False)
m['sensitivity_direction_agrees'] = (m.log2FC_primary * m.log2FC_sensitivity > 0).fillna(False)
m['exclude_index_gene'] = m.index.to_series().str.split('.').str[0].eq('ENSG00000198911').values
m['annotation_pass'] = m.annotation_status.eq('exact_match')
m['signature_selected'] = (m.primary_fdr_pass & m.primary_effect_pass & m.fbs_direction_agrees &
                           m.sensitivity_direction_agrees & ~m.exclude_index_gene & m.annotation_pass)
m['direction'] = 'not_selected'
m.loc[m.signature_selected & (m.log2FC_primary < 0), 'direction'] = 'decreased_in_KO'
m.loc[m.signature_selected & (m.log2FC_primary > 0), 'direction'] = 'increased_in_KO'
m['gene_id'] = m.index
cols = ['gene_id', 'gene_symbol', 'gene_type', 'annotation_status', 'baseMean', 'log2FC_primary',
        'padj_primary', 'log2FC_fbs', 'padj_fbs', 'log2FC_sensitivity', 'padj_sensitivity',
        'primary_fdr_pass', 'primary_effect_pass', 'fbs_direction_agrees',
        'sensitivity_direction_agrees', 'exclude_index_gene', 'annotation_pass',
        'signature_selected', 'direction']
m = m[cols]
selected = m[m.signature_selected].reset_index(drop=True).sort_values(['padj_primary', 'gene_id'], ignore_index=True)
positive = selected[selected.direction == 'decreased_in_KO']
negative = selected[selected.direction == 'increased_in_KO']
if selected.empty:
    raise ValueError('No genes selected; inspect filters before proceeding')
OUT.mkdir(parents=True, exist_ok=False)
outputs = {
    'SIGNATURE_SELECTION_AUDIT.tsv': m,
    'SREBF2_SIGNATURE_ALL.tsv': selected,
    'SREBF2_SIGNATURE_POSITIVE.tsv': positive,
    'SREBF2_SIGNATURE_NEGATIVE.tsv': negative,
}
for filename, frame in outputs.items():
    frame.to_csv(OUT / filename, sep='\t', index=False, float_format='%.12g')
manifest = {
    'signature_version': 'discovery_v1', 'date_created': str(date.today()),
    'classification': 'perturbation-derived candidate signature; binding not established',
    'source_dataset': 'GSE267018', 'validation_data_used_for_selection': False,
    'reference_assembly': 'GRCh38', 'gene_annotation': 'GENCODE v47',
    'rules': ['primary: KO vs WT under depletion, padj < 0.05, abs(log2FC) >= 1',
              'same nonzero effect direction in FBS and SRR28966285-excluded sensitivity contrasts',
              'exact-match annotation only', 'exclude SREBF2 itself',
              'no FBS/sensitivity significance requirement'],
    'scoring_policy': 'NOT FROZEN; do not evaluate validation until scoring, mapping and missing-gene rules are specified',
    'counts': {'genes_evaluated': len(m), 'selected': len(selected),
               'decreased_in_KO': len(positive), 'increased_in_KO': len(negative)},
    'inputs_sha256': {name: sha256(path) for name, path in SOURCES.items()},
    'outputs_sha256': {name: sha256(OUT / name) for name in outputs},
}
(OUT / 'SIGNATURE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
(OUT / 'METHODS.md').write_text('''# SREBF2 candidate signature — discovery v1\n\nDiscovery only (GSE267018); no GSE271000 data were read.\n\nSelection: primary KO-vs-WT cholesterol-depletion BH-adjusted p < 0.05, |log2FC| >= 1, concordant nonzero direction under FBS and in the exclusion-of-SRR28966285 sensitivity fit; exact-match gene annotation; exclude SREBF2. FBS and sensitivity significance are **not** required.\n\nA negative KO-vs-WT fold change means expression decreased following SREBF2 knockout (candidate positively associated gene). A positive KO-vs-WT fold change means increased expression following knockout (candidate negatively associated gene). These are not proven direct binding targets.\n\nGene membership and direction are frozen here. A fully specified scoring function, ID mapping, coverage and missing-gene policy must be frozen before validation. Do not use validation DE to alter this gene list.\n''')
print(json.dumps(manifest['counts'], indent=2))
print('Frozen candidate gene membership at:', OUT)
