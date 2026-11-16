#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase12/phase12_expression_stratified.py
# Purpose     : Phase 12K: expression-decile and nonlinear-expression sensitivity
#               analyses
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 12K: expression-decile and nonlinear-expression sensitivity analyses.

Run from repository root in the same environment as Phase 12J.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

SOURCE = Path('results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv')
OUT = Path('results/phase12')
OUT.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(SOURCE, sep='\t', low_memory=False)
needed = ['_merge','signature_selected','direction','baseMean','gene_length','gene_type',
          'mean_promoter_H3K27ac','mean_promoter_SREBP2_CholDep','mean_promoter_Input_CholDep']
missing = sorted(set(needed) - set(df.columns))
if missing:
    raise ValueError(f'Missing columns: {missing}')
df = df.loc[df['_merge'].eq('both')].copy()
df['signature'] = df['signature_selected'].astype(str).str.lower().eq('true')
df['baseMean'] = pd.to_numeric(df['baseMean'], errors='coerce')
df['gene_length'] = pd.to_numeric(df['gene_length'], errors='coerce')
if (df['baseMean'].dropna() < 0).any() or (df['gene_length'].dropna() <= 0).any():
    raise ValueError('Invalid expression or gene length')
df['log_baseMean'] = np.log10(df['baseMean'] + 1)
df['log_gene_length'] = np.log10(df['gene_length'])
df['biotype_group'] = df['gene_type'].where(df['gene_type'].isin(['protein_coding','lncRNA']), 'other')
for key, col in {'H3K27ac':'mean_promoter_H3K27ac',
                 'SREBP2_CholDep':'mean_promoter_SREBP2_CholDep',
                 'Input_CholDep':'mean_promoter_Input_CholDep'}.items():
    v = pd.to_numeric(df[col], errors='coerce')
    if (v.dropna() < 0).any():
        raise ValueError(f'Negative signal in {col}')
    v = np.log1p(v)
    df['z_' + key] = (v - v.mean()) / v.std()
required = ['log_baseMean','log_gene_length','z_H3K27ac','z_SREBP2_CholDep','z_Input_CholDep']
df = df.replace([np.inf,-np.inf], np.nan).dropna(subset=required).copy()
df['expression_decile'] = pd.qcut(df['log_baseMean'], q=10, labels=False, duplicates='drop').astype(int) + 1
if df['expression_decile'].nunique() != 10:
    raise RuntimeError('Could not construct ten expression deciles')
print('========== PHASE 12K ==========')
print('Matched genes:',len(df),'Signature genes:',int(df['signature'].sum()))
print('Expression strata:',df['expression_decile'].nunique())
comparisons = {
    'decreased_vs_background':df.loc[(~df['signature']) | df['direction'].eq('decreased_in_KO')].copy(),
    'increased_vs_background':df.loc[(~df['signature']) | df['direction'].eq('increased_in_KO')].copy(),
    'increased_vs_decreased':df.loc[df['signature']].copy()
}
results, counts = [], []
for comparison, sub in comparisons.items():
    if comparison == 'increased_vs_decreased':
        sub['outcome'] = sub['direction'].eq('increased_in_KO').astype(int)
    else:
        sub['outcome'] = sub['signature'].astype(int)
    for population, d in [('all_biotypes',sub),('protein_coding',sub.loc[sub['gene_type'].eq('protein_coding')].copy())]:
        for decile, group in d.groupby('expression_decile'):
            counts.append(dict(comparison=comparison,population=population,expression_decile=int(decile),
                               n=len(group),cases=int(group['outcome'].sum()),
                               controls=int((1-group['outcome']).sum()),
                               median_baseMean=float(group['baseMean'].median())))
        base = 'z_H3K27ac + z_SREBP2_CholDep + z_Input_CholDep + log_gene_length'
        if population == 'all_biotypes':
            base += ' + C(biotype_group)'
        specs = {
            'expression_deciles':'outcome ~ ' + base + ' + C(expression_decile)',
            'expression_spline':'outcome ~ ' + base + ' + bs(log_baseMean, df=4, degree=3)'
        }
        for spec, formula in specs.items():
            if spec == 'expression_deciles':
                grouped = d.groupby('expression_decile')['outcome'].agg(['sum','count'])
                sparse = grouped.loc[(grouped['sum']==0) | (grouped['sum']==grouped['count'])]
                if len(sparse):
                    print(f'WARNING: {comparison}/{population} has single-outcome deciles: {sparse.index.tolist()}')
            try:
                fit = smf.logit(formula, data=d).fit(disp=False,maxiter=300)
            except Exception as exc:
                print(f'ERROR fitting {comparison}/{population}/{spec}: {exc}')
                raise
            converged = bool(fit.mle_retvals['converged'])
            print(f'{comparison} | {population} | {spec} | n={int(fit.nobs)} | converged={converged}')
            ci = fit.conf_int()
            for predictor in ['z_H3K27ac','z_SREBP2_CholDep']:
                results.append(dict(comparison=comparison,population=population,specification=spec,
                                    predictor=predictor,n=int(fit.nobs),cases=int(d['outcome'].sum()),
                                    odds_ratio=float(np.exp(fit.params[predictor])),
                                    ci_lower=float(np.exp(ci.loc[predictor,0])),
                                    ci_upper=float(np.exp(ci.loc[predictor,1])),
                                    p_value=float(fit.pvalues[predictor]),converged=converged))
models = pd.DataFrame(results)
models['fdr_bh'] = multipletests(models['p_value'],method='fdr_bh')[1]
strata = pd.DataFrame(counts)
models.to_csv(OUT/'PHASE12_EXPRESSION_STRATIFIED_MODELS.tsv',sep='\t',index=False)
strata.to_csv(OUT/'PHASE12_EXPRESSION_STRATA_COUNTS.tsv',sep='\t',index=False)
print('\n========== MODEL RESULTS ==========')
print(models.to_string(index=False))
print('\n========== STRATA COUNTS ==========')
print(strata.to_string(index=False))
print('\nPHASE 12K COMPLETE')
