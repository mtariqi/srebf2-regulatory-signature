#!/usr/bin/env python3
"""Frozen SREBF2 signature concordance in GEO GSE271000 processed DESeq2 contrasts."""
from pathlib import Path
import gzip, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path.cwd()
SIG = ROOT / 'results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv'
RAW = ROOT / 'data/raw/validation/GSE271000'
OUT = ROOT / 'results/validation/GSE271000/processed_DE_v1'

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

if OUT.exists():
    raise SystemExit(f'Output already exists; refusing to overwrite: {OUT}')
sig = pd.read_csv(SIG, sep='\t', dtype={'gene_id': str})
required = {'gene_id','gene_symbol','log2FC_primary','signature_selected','direction'}
if not required.issubset(sig.columns):
    raise SystemExit(f'Signature missing columns: {required - set(sig.columns)}')
if len(sig) != 1326 or not sig.gene_id.is_unique:
    raise SystemExit('Frozen signature count or gene_id uniqueness invariant failed')
if not sig.signature_selected.astype(str).str.lower().eq('true').all():
    raise SystemExit('Signature contains unselected genes')
sig['expected_sign'] = np.sign(pd.to_numeric(sig.log2FC_primary, errors='coerce'))
if sig.expected_sign.isna().any() or sig.expected_sign.eq(0).any():
    raise SystemExit('Invalid discovery log2FC')
if not ((sig.direction.eq('decreased_in_KO') & sig.expected_sign.eq(-1)) |
        (sig.direction.eq('increased_in_KO') & sig.expected_sign.eq(1))).all():
    raise SystemExit('Direction inconsistent with discovery log2FC')
symbol_counts = sig.gene_symbol.value_counts(dropna=False)
sig['unique_signature_symbol'] = (sig.gene_symbol.notna() & sig.gene_symbol.astype(str).ne('') &
                                   sig.gene_symbol.map(symbol_counts).eq(1))

inputs = {str(SIG): sha256(SIG)}
rows = []
summaries = []
for hairpin in ('664','665','667'):
    path = RAW / f'GSE271000_DESeq_HepG2_shSREBP2_{hairpin}.txt.gz'
    inputs[str(path)] = sha256(path)
    with gzip.open(path, 'rt') as handle:
        de = pd.read_csv(handle, sep='\t', index_col=0)
    de.index = de.index.astype(str)
    needed = {'baseMean','log2FoldChange','pvalue','padj'}
    if not needed.issubset(de.columns):
        raise SystemExit(f'{path.name}: missing {needed - set(de.columns)}')
    duplicates = set(de.index[de.index.duplicated(keep=False)])
    de = de.loc[~de.index.isin(duplicates), list(needed)].copy()
    de = de.rename(columns={'baseMean':'validation_baseMean', 'pvalue':'validation_pvalue',
                            'padj':'validation_padj'})
    mapped = sig.copy()
    mapped['hairpin'] = hairpin
    mapped['mapping_eligible'] = mapped.unique_signature_symbol & ~mapped.gene_symbol.isin(duplicates)
    mapped = mapped.join(de, on='gene_symbol')
    value_cols = ['validation_baseMean','log2FoldChange','validation_pvalue','validation_padj']
    mapped.loc[~mapped.mapping_eligible, value_cols] = np.nan
    mapped['validation_sign'] = np.sign(pd.to_numeric(mapped.log2FoldChange, errors='coerce'))
    evaluable = mapped.mapping_eligible & mapped.validation_sign.notna() & mapped.validation_sign.ne(0)
    mapped['concordant'] = pd.Series(pd.NA, index=mapped.index, dtype='boolean')
    mapped.loc[evaluable, 'concordant'] = mapped.loc[evaluable, 'expected_sign'].eq(mapped.loc[evaluable, 'validation_sign']).to_numpy()
    valid = mapped.loc[evaluable]
    n = len(valid)
    k = int(valid.concordant.sum())
    rho = float(spearmanr(valid.log2FC_primary, valid.log2FoldChange).statistic) if n >= 3 else float('nan')
    summaries.append({'hairpin':hairpin,'signature_total':len(sig),
                      'mapped_with_nonzero_LFC':n,'concordant':k,
                      'concordance_fraction':k/n if n else np.nan,
                      'spearman_rho':rho,
                      'binomial_p_exploratory':binomtest(k,n,0.5,alternative='greater').pvalue if n else np.nan,
                      'ambiguous_signature_symbols':int((~sig.unique_signature_symbol).sum()),
                      'duplicate_validation_symbols':len(duplicates)})
    rows.append(mapped)
full = pd.concat(rows, ignore_index=True)
summary = pd.DataFrame(summaries)
wide = full.pivot(index='gene_id', columns='hairpin', values='concordant').reindex(columns=['664','665','667'])
wide['all_three_evaluable'] = wide[['664','665','667']].notna().all(axis=1)
wide['all_three_concordant'] = wide.all_three_evaluable & wide[['664','665','667']].eq(True).all(axis=1)
wide['concordant_count'] = wide[['664','665','667']].eq(True).sum(axis=1)
OUT.mkdir(parents=True)
full.to_csv(OUT/'VALIDATION_GENE_MAPPING.tsv',sep='\t',index=False)
summary.to_csv(OUT/'VALIDATION_DIRECTION_CONCORDANCE.tsv',sep='\t',index=False)
wide.reset_index().to_csv(OUT/'VALIDATION_THREE_SHRNA_CONSENSUS.tsv',sep='\t',index=False)
manifest = {'analysis':'GSE271000_processed_DE_v1','input_sha256':inputs,
            'frozen_signature_size':1326,'assumed_contrast':'shSREBP2_vs_shControl',
            'contrast_orientation_verified':False,
            'mapping':'exact unique gene symbols; ambiguous excluded',
            'limitations':['No sample-level scores','Contrast orientation not independently verified',
                           'Gene dependence makes binomial p values exploratory',
                           'Cell type and perturbation differ from discovery']}
(OUT/'VALIDATION_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Per-hairpin concordance:')
print(summary.to_string(index=False))
print('\nThree-shRNA consensus:')
print('All three evaluable:',int(wide.all_three_evaluable.sum()))
print('All three concordant:',int(wide.all_three_concordant.sum()))
print('\nSaved:',OUT)
