#!/usr/bin/env python3
"""Exploratory adjusted promoter analysis, preserves all frozen v1 files.
Run from repository root: conda run -n srebf2-binding python scripts/srebf2_promoter_adjusted_v2.py
Requires pandas, numpy, scipy, statsmodels.
"""
import csv, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf

BASE=Path('results/binding/GSE282800')
SRC=BASE/'SREBF2_DISCOVERY_UNIVERSE_PROMOTER_SIGNAL_v1.tsv'
OUT=BASE/'adjusted_v2'
OUT.mkdir(parents=True, exist_ok=True)
SEED=20261011
rng=np.random.default_rng(SEED)

raw=pd.read_csv(SRC,sep='\t',dtype={'gene_id':str},low_memory=False)
if len(raw)!=17158 or raw.gene_id.duplicated().any():
    raise ValueError('Frozen universe count or unique gene IDs failed')
needed=['mapping_status','primary_nuclear_chromosome','signature_selected','direction',
        'gene_type','baseMean','FBS_mean','depletion_mean','depletion_minus_FBS','chrom']
missing=set(needed)-set(raw.columns)
if missing: raise ValueError(f'Missing columns: {sorted(missing)}')

eligible=raw.loc[(raw.mapping_status=='mapped_with_signal') &
                 (raw.primary_nuclear_chromosome.astype(str).str.lower()=='true')].copy()
eligible['group']=np.where(eligible.signature_selected.astype(str).str.lower()=='true',
                           eligible.direction,'unselected')
expected={'decreased_in_KO':789,'increased_in_KO':527,'unselected':15781}
if eligible.group.value_counts().to_dict()!=expected:
    raise ValueError(f'Unexpected v1 groups: {eligible.group.value_counts().to_dict()}')
for field in ['baseMean','FBS_mean','depletion_mean','depletion_minus_FBS']:
    eligible[field]=pd.to_numeric(eligible[field],errors='coerce')
eligible=eligible.loc[np.isfinite(eligible.baseMean) & (eligible.baseMean>=0) &
                      np.isfinite(eligible.depletion_minus_FBS) &
                      eligible.gene_type.notna()].copy()
eligible['log_baseMean']=np.log1p(eligible.baseMean)
# Group uncommon biotypes into an explicit category, to avoid unstable dummy coefficients.
counts=eligible.gene_type.value_counts()
eligible['biotype_group']=eligible.gene_type.where(eligible.gene_type.map(counts)>=100,'other_rare')
eligible['group']=pd.Categorical(eligible.group,
    categories=['unselected','decreased_in_KO','increased_in_KO'])

# HC3 protects against heteroskedasticity, NOT spatial correlation among nearby genes.
formula=('depletion_minus_FBS ~ C(group, Treatment(reference="unselected")) '
         '+ log_baseMean + I(log_baseMean ** 2) + C(biotype_group)')
model=smf.ols(formula,data=eligible).fit(cov_type='HC3')
reg=[]
for term in model.params.index:
    if 'C(group' in term:
        ci=model.conf_int().loc[term]
        reg.append({'term':term,'coefficient':model.params[term],
                    'HC3_CI_low':ci.iloc[0],'HC3_CI_high':ci.iloc[1],
                    'p_exploratory':model.pvalues[term]})
pd.DataFrame(reg).to_csv(OUT/'ADJUSTED_REGRESSION_v2.tsv',sep='\t',index=False)
(OUT/'ADJUSTED_REGRESSION_SUMMARY_v2.txt').write_text(model.summary().as_text()+'\n')

# Matched comparison: each selected gene compared with ONE unselected gene in the
# same biotype and within the same global log-expression decile. Sampling without
# replacement per group; separate negative/positive comparisons can reuse controls.
eligible['expression_bin']=pd.qcut(eligible.log_baseMean.rank(method='first'),10,labels=False)
matched=[]
coverage=[]
for target in ['decreased_in_KO','increased_in_KO']:
    cases=eligible.loc[eligible.group==target].copy()
    controls=eligible.loc[eligible.group=='unselected'].copy()
    selected_case=[]; selected_control=[]
    for (bio,bin_id),part in cases.groupby(['biotype_group','expression_bin'],observed=True):
        pool=controls.loc[(controls.biotype_group==bio)&(controls.expression_bin==bin_id)]
        k=min(len(part),len(pool))
        if k:
            # Fixed seed, deterministic sampling within each stratum.
            selected_case.extend(rng.choice(part.index.to_numpy(),size=k,replace=False).tolist())
            selected_control.extend(rng.choice(pool.index.to_numpy(),size=k,replace=False).tolist())
    x=eligible.loc[selected_case,'depletion_minus_FBS'].to_numpy()
    y=eligible.loc[selected_control,'depletion_minus_FBS'].to_numpy()
    if len(x)<2: raise ValueError(f'Insufficient matched genes: {target}')
    u,p=mannwhitneyu(x,y,alternative='two-sided',method='asymptotic')
    matched.append({'group':target,'matched_cases':len(x),'total_cases':len(cases),
        'matched_controls':len(y),'median_case':np.median(x),'median_control':np.median(y),
        'median_difference':np.median(x)-np.median(y),'AUC':u/(len(x)*len(y)),
        'p_exploratory_unpaired':p})
    coverage.append({'group':target,'case_coverage_fraction':len(x)/len(cases),
                     'unmatched_cases':len(cases)-len(x)})
    pd.DataFrame({'case_gene_id':eligible.loc[selected_case,'gene_id'].to_numpy(),
                  'control_gene_id':eligible.loc[selected_control,'gene_id'].to_numpy()}).to_csv(
                      OUT/f'MATCHED_PAIRS_{target}_v2.tsv',sep='\t',index=False)
pd.DataFrame(matched).to_csv(OUT/'EXPRESSION_BIOTYPE_MATCHED_v2.tsv',sep='\t',index=False)

# Protein-coding-only sensitivity; keep all 3 groups and provide effect sizes.
protein=eligible.loc[eligible.gene_type=='protein_coding']
protein_results=[]
for target in ['decreased_in_KO','increased_in_KO']:
    x=protein.loc[protein.group==target,'depletion_minus_FBS'].to_numpy()
    y=protein.loc[protein.group=='unselected','depletion_minus_FBS'].to_numpy()
    if len(x) and len(y):
        u,p=mannwhitneyu(x,y,alternative='two-sided',method='asymptotic')
        protein_results.append({'group':target,'n_case':len(x),'n_unselected':len(y),
            'median_case':np.median(x),'median_unselected':np.median(y),
            'AUC':u/(len(x)*len(y)),'p_exploratory':p})
pd.DataFrame(protein_results).to_csv(OUT/'PROTEIN_CODING_ONLY_v2.tsv',sep='\t',index=False)

sha=hashlib.sha256(SRC.read_bytes()).hexdigest()
manifest={'analysis':'adjusted_v2','source':str(SRC),'source_sha256':sha,
    'seed':SEED,'input_rows':len(raw),'eligible_before_covariate_filter':17097,
    'eligible_after_covariate_filter':len(eligible),
    'excluded_due_to_missing_covariates':17097-len(eligible),
    'gene_type_counts':eligible.gene_type.value_counts().to_dict(),
    'model':formula,'covariance':'HC3; does not correct genomic dependence',
    'matching':'same grouped biotype and global expression decile; 1:1 without replacement within comparison',
    'matching_coverage':coverage,
    'caveats':['Exploratory, noncausal, no replicate-level differential binding',
               'Single random matched draw; matching p-values are descriptive only',
               'No genomic-block uncertainty estimation',
               'Gene-level TSS and promoter intervals may overlap']}
(OUT/'ADJUSTED_ANALYSIS_MANIFEST_v2.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Source SHA256:',sha)
print('Eligible genes:',len(eligible),'of 17097; dropped for covariates:',17097-len(eligible))
print('\nAdjusted regression (HC3, exploratory):')
print(pd.DataFrame(reg).to_string(index=False))
print('\nMatched comparisons (single fixed-seed draw):')
print(pd.DataFrame(matched).to_string(index=False))
print('\nProtein-coding-only sensitivity:')
print(pd.DataFrame(protein_results).to_string(index=False))
print('\nSaved:',OUT)
