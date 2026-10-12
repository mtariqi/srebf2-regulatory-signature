#!/usr/bin/env python3
"""Sensitivity audit of SREBF2 promoter-signal association; never modifies v1/v2."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

BASE=Path('results/binding/GSE282800')
SRC=BASE/'SREBF2_DISCOVERY_UNIVERSE_PROMOTER_SIGNAL_v1.tsv'
OUT=BASE/'robustness_v3'; OUT.mkdir(parents=True,exist_ok=True)
EXPECTED='a3cf94930db86e82ddf18e58f59aa870d92d0c94355cf8e29eaecf2891e4d947'
sha=hashlib.sha256(SRC.read_bytes()).hexdigest()
if sha!=EXPECTED: raise ValueError(f'Unexpected source SHA256: {sha}')
df=pd.read_csv(SRC,sep='\t',low_memory=False)
if len(df)!=17158 or df.gene_id.duplicated().any(): raise ValueError('Universe integrity failure')
df=df[(df.mapping_status=='mapped_with_signal') &
      (df.primary_nuclear_chromosome.astype(str).str.lower()=='true')].copy()
df['group']=np.where(df.signature_selected.astype(str).str.lower()=='true',df.direction,'unselected')
if df.group.value_counts().to_dict()!={'unselected':15781,'decreased_in_KO':789,'increased_in_KO':527}:
    raise ValueError('Unexpected group counts')
for c in ['baseMean','depletion_minus_FBS']:
    df[c]=pd.to_numeric(df[c],errors='coerce')
df=df[np.isfinite(df.baseMean)&(df.baseMean>=0)&np.isfinite(df.depletion_minus_FBS)&df.gene_type.notna()].copy()
df['log_baseMean']=np.log1p(df.baseMean)
counts=df.gene_type.value_counts()
df['biotype_group']=df.gene_type.where(df.gene_type.map(counts)>=100,'other_rare')
df['group']=pd.Categorical(df.group,categories=['unselected','decreased_in_KO','increased_in_KO'])
formula='OUTCOME ~ C(group, Treatment(reference="unselected")) + log_baseMean + I(log_baseMean ** 2) + C(biotype_group)'
y=df.depletion_minus_FBS
lo,hi=np.quantile(y,[0.01,0.99])
lo5,hi5=np.quantile(y,[0.05,0.95])
scale=float(np.median(np.abs(y-np.median(y))))
if scale<=0: scale=float(np.std(y))
variants={
 'raw':y,
 'winsor_1pct':y.clip(lo,hi),
 'winsor_5pct':y.clip(lo5,hi5),
 'asinh_MAD_scaled':np.arcsinh(y/scale),
}
records=[]
for label,values in variants.items():
    temp=df.copy();temp['OUTCOME']=values
    fit=smf.ols(formula,data=temp).fit(cov_type='HC3')
    for term in fit.params.index:
        if 'C(group' not in term:continue
        ci=fit.conf_int().loc[term]
        records.append({'model':label,'term':term,'coefficient':fit.params[term],
                        'CI_low_HC3':ci.iloc[0],'CI_high_HC3':ci.iloc[1],
                        'p_exploratory':fit.pvalues[term],'R_squared':fit.rsquared,
                        'outcome_scale':('asinh(delta/MAD)' if label=='asinh_MAD_scaled' else 'bigWig delta')})
res=pd.DataFrame(records)
res.to_csv(OUT/'ROBUST_REGRESSION_SENSITIVITY_v3.tsv',sep='\t',index=False)
# Outlier concentration: extreme genes, group shares and their contribution to group means.
rows=[]
for group,part in df.groupby('group',observed=True):
    v=part.depletion_minus_FBS.to_numpy()
    top=np.sort(v)[::-1]
    rows.append({'group':group,'n':len(v),'median':np.median(v),'mean':np.mean(v),
                 'p95':np.quantile(v,.95),'p99':np.quantile(v,.99),
                 'max':np.max(v),'min':np.min(v),
                 'mean_top_1pct_removed':np.mean(np.sort(v)[:max(1,int(np.floor(.99*len(v))))])})
pd.DataFrame(rows).to_csv(OUT/'OUTLIER_DISTRIBUTIONS_v3.tsv',sep='\t',index=False)
# The v2 case-control file pairs genes arbitrarily within strata; do not treat them as matched by distance.
# Summarize repeated matching using 100 independent random control draws within the same biotype/bin.
df['expression_bin']=pd.qcut(df.log_baseMean.rank(method='first'),10,labels=False)
rng=np.random.default_rng(20261011)
match_rows=[]
for target in ['decreased_in_KO','increased_in_KO']:
    cases=df[df.group==target]
    controls=df[df.group=='unselected']
    pools=[]
    for key,part in cases.groupby(['biotype_group','expression_bin'],observed=True):
        pool=controls[(controls.biotype_group==key[0])&(controls.expression_bin==key[1])]
        k=min(len(part),len(pool))
        if k:pools.append((part.depletion_minus_FBS.to_numpy(),pool.depletion_minus_FBS.to_numpy(),k))
    for draw in range(100):
        casevals=[];controlvals=[]
        for cv,pv,k in pools:
            casevals.extend(rng.choice(cv,k,replace=False))
            controlvals.extend(rng.choice(pv,k,replace=False))
        casevals=np.asarray(casevals);controlvals=np.asarray(controlvals)
        match_rows.append({'group':target,'draw':draw+1,'n':len(casevals),
            'difference_of_medians':np.median(casevals)-np.median(controlvals),
            'difference_of_means':np.mean(casevals)-np.mean(controlvals)})
match=pd.DataFrame(match_rows)
match.to_csv(OUT/'REPEATED_MATCHING_100_DRAWS_v3.tsv',sep='\t',index=False)
print('Regression sensitivity (HC3 exploratory):')
print(res[['model','term','coefficient','CI_low_HC3','CI_high_HC3']].to_string(index=False))
print('\nGroup distribution summary:')
print(pd.DataFrame(rows).to_string(index=False))
print('\nRepeated matching: range across 100 draws (NOT confidence intervals):')
print(match.groupby('group')['difference_of_medians'].agg(['min','median','max']).to_string())
manifest={'source':str(SRC),'source_sha256':sha,'analysis':'robustness_v3',
 'n_eligible':len(df),'winsor_bounds_1pct':[float(lo),float(hi)],
 'winsor_bounds_5pct':[float(lo5),float(hi5)],'asinh_MAD_scale':scale,
 'matching_draws':100,'seed':20261011,
 'limitations':['Exploratory association, not differential binding',
 'HC3 does not account for promoter overlap or genomic dependence',
 'Repeated matching range reflects random draws, not sampling uncertainty',
 'Outcome transformations change coefficient scales',
 'Expression deciles and grouped biotypes do not ensure complete confounder balance']}
(OUT/'ROBUSTNESS_MANIFEST_v3.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('\nSaved:',OUT)
