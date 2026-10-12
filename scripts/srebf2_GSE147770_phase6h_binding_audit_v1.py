#!/usr/bin/env python3
"""Phase 6H: exploratory cross-study promoter-signal audit; run at repository root.
Requires pandas numpy scipy; reads frozen Phase 6G and GSE282800 outputs.
GSE282800 LPPS is cholesterol depletion, NOT HAP1 -FBS serum starvation.
No ChIP peak calls or direct-binding claims are made.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT=Path.cwd()
P_BIND=Path('results/binding/GSE282800/SREBF2_DISCOVERY_UNIVERSE_PROMOTER_SIGNAL_v1.tsv')
P_LEAD=Path('results/validation/GSE147770/mitochondrial_robustness_v1/OXPHOS_13_GENES_v1.tsv')
P_AUDIT=Path('results/validation/GSE147770/direction_audit_v1/SIGNATURE_DIRECTION_AUDIT_v1.tsv')
P_WALD=Path('results/validation/GSE147770/interaction_v1/interaction_Wald_v1.tsv')
OUT=Path('results/validation/GSE147770/binding_support_v1')
for p in [P_BIND,P_LEAD,P_AUDIT,P_WALD]:
    if not p.is_file():raise SystemExit(f'Missing input: {p}')
if OUT.exists():raise SystemExit(f'Refusing to overwrite existing directory: {OUT}')

def load(path):return pd.read_csv(path,sep='\t',low_memory=False)
def unversioned(s):return s.astype(str).str.replace(r'\.\d+$','',regex=True)
def bool_col(s):return s.astype(str).str.lower().isin(['true','1','t'])
def desc(x):
    v=pd.to_numeric(x,errors='coerce').dropna().to_numpy(float)
    return {'n':len(v),'median':float(np.median(v)) if len(v) else None,'mean':float(np.mean(v)) if len(v) else None}
def compare(case,ctrl):
    a=pd.to_numeric(case,errors='coerce').dropna().to_numpy(float)
    b=pd.to_numeric(ctrl,errors='coerce').dropna().to_numpy(float)
    if not len(a) or not len(b):return dict(n_case=len(a),n_control=len(b),median_case=np.nan,median_control=np.nan,median_difference=np.nan,auc=np.nan,p_two_sided=np.nan)
    u,p=mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
    return dict(n_case=len(a),n_control=len(b),median_case=np.median(a),median_control=np.median(b),median_difference=np.median(a)-np.median(b),auc=u/(len(a)*len(b)),p_two_sided=p)

bind=load(P_BIND);lead=load(P_LEAD);audit=load(P_AUDIT);wald=load(P_WALD)
for name,df,required in [
 ('binding',bind,['gene_id','mapping_status','primary_nuclear_chromosome','FBS_mean','depletion_mean','depletion_minus_FBS','baseMean','gene_type']),
 ('leading',lead,['unversioned_id','plus_log2FC','minus_log2FC','interaction_padj']),
 ('audit',audit,['unversioned_id','direction','interaction_FDR05']),
 ('wald',wald,['unversioned_id','padj'])]:
    missing=set(required)-set(df.columns)
    if missing:raise SystemExit(f'{name} missing columns: {sorted(missing)}')
if len(lead)!=13 or lead.unversioned_id.duplicated().any():raise SystemExit('Expected exactly 13 unique Phase 6G leading genes')
if len(bind)!=17158 or bind.gene_id.duplicated().any():raise SystemExit('Unexpected binding universe rows or duplicate gene_id')
bind['unversioned_id']=unversioned(bind.gene_id)
if bind.unversioned_id.duplicated().any():raise SystemExit('Version-stripped binding gene IDs are not unique')
for c in ['FBS_mean','depletion_mean','depletion_minus_FBS','baseMean']:
    bind[c]=pd.to_numeric(bind[c],errors='coerce')
# Binding eligible universe is GSE267018 discovery gene universe with usable promoter coverage.
bind['promoter_eligible']=(bind.mapping_status.eq('mapped_with_signal') & bool_col(bind.primary_nuclear_chromosome) & np.isfinite(bind.FBS_mean) & np.isfinite(bind.depletion_mean) & np.isfinite(bind.depletion_minus_FBS))
wald_ids=set(wald.loc[wald.padj.notna(),'unversioned_id'].astype(str))
audit_ids=set(audit.loc[audit.unversioned_id.isin(wald_ids),'unversioned_id'].astype(str))
if len(wald_ids)!=14171 or len(audit_ids)!=807:raise SystemExit(f'Unexpected HAP1 evaluable counts: {len(wald_ids)} / signature {len(audit_ids)}')
if not set(lead.unversioned_id).issubset(audit_ids):raise SystemExit('Leading genes not all in evaluable signature')
joined=lead.merge(bind,how='left',on='unversioned_id',validate='one_to_one',suffixes=('','_discovery'))
joined['promoter_eligible']=joined.promoter_eligible.fillna(False)
joined['has_any_promoter_signal']=joined.promoter_eligible & ((joined.FBS_mean>0)|(joined.depletion_mean>0))
# Frozen-signature background: genes with HAP1 interaction padj AND promoter coverage in GSE282800.
background=bind[bind.unversioned_id.isin(audit_ids) & bind.promoter_eligible].copy()
case=joined[joined.promoter_eligible].copy()
control=background[~background.unversioned_id.isin(lead.unversioned_id)].copy()
if len(case)<2 or len(control)<10:raise SystemExit('Insufficient overlap with promoter measurements')
OUT.mkdir(parents=True)
joined.to_csv(OUT/'OXPHOS_13_PROMOTER_SIGNAL_v1.tsv',sep='\t',index=False)
rows=[]
for metric in ['FBS_mean','depletion_mean','depletion_minus_FBS']:
    rows.append(dict(metric=metric,background='other_evaluable_signature',**compare(case[metric],control[metric])))
# Secondary comparison to other interaction-significant signature genes.
interaction_ids=set(audit.loc[bool_col(audit.interaction_FDR05)&audit.unversioned_id.isin(wald_ids),'unversioned_id'])
control_int=control[control.unversioned_id.isin(interaction_ids)]
for metric in ['FBS_mean','depletion_mean','depletion_minus_FBS']:
    rows.append(dict(metric=metric,background='other_interaction_significant_signature',**compare(case[metric],control_int[metric])))
pd.DataFrame(rows).to_csv(OUT/'PROMOTER_SIGNAL_COMPARISONS_v1.tsv',sep='\t',index=False)
# Expression and biotype matched sampling. Cases with no available stratum controls are excluded explicitly.
combined=pd.concat([case.assign(group='case'),control.assign(group='control')],ignore_index=True)
combined=combined[np.isfinite(combined.baseMean)&(combined.baseMean>=0)&combined.gene_type.notna()].copy()
combined['expr_bin']=pd.qcut(np.log1p(combined.baseMean).rank(method='first'),q=5,labels=False,duplicates='drop')
combined['stratum']=combined.gene_type.astype(str)+'|'+combined.expr_bin.astype(str)
case_m=combined[combined.group=='case'];ctrl_m=combined[combined.group=='control']
matched_case=case_m[case_m.stratum.isin(ctrl_m.stratum)].copy()
matched_controls={s:part for s,part in ctrl_m.groupby('stratum')}
rng=np.random.default_rng(20261012)
match_rows=[]
for draw in range(500):
    sampled=[]
    for s,part in matched_case.groupby('stratum'):
        pool=matched_controls[s]
        sampled.append(pool.iloc[rng.choice(len(pool),size=len(part),replace=len(pool)<len(part))])
    matched=pd.concat(sampled,ignore_index=True) if sampled else ctrl_m.iloc[0:0]
    for metric in ['FBS_mean','depletion_mean','depletion_minus_FBS']:
        x=matched_case[metric].to_numpy(float);y=matched[metric].to_numpy(float)
        match_rows.append({'draw':draw+1,'metric':metric,'n_matched_cases':len(x),'n_sampled_controls':len(y),'median_difference':float(np.median(x)-np.median(y)) if len(x) else np.nan,'mean_difference':float(np.mean(x)-np.mean(y)) if len(x) else np.nan})
match=pd.DataFrame(match_rows)
match.to_csv(OUT/'MATCHED_SENSITIVITY_500_DRAWS_v1.tsv',sep='\t',index=False)
summary=match.groupby('metric',as_index=False).agg(n_matched_cases=('n_matched_cases','first'),median_of_median_differences=('median_difference','median'),min_median_difference=('median_difference','min'),max_median_difference=('median_difference','max'))
summary.to_csv(OUT/'MATCHED_SENSITIVITY_SUMMARY_v1.tsv',sep='\t',index=False)
# Nonzero promoter signal is descriptive, not a peak or binding threshold.
coverage=pd.DataFrame([{'leading_genes':len(lead),'leading_with_promoter_coverage':len(case),'leading_nonzero_either_track':int(joined.has_any_promoter_signal.sum()),'evaluable_signature_with_promoter_coverage':len(background),'comparison_signature_genes':len(control),'comparison_interaction_significant_genes':len(control_int),'matched_leading_genes':len(matched_case)}])
coverage.to_csv(OUT/'PROMOTER_COVERAGE_v1.tsv',sep='\t',index=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'phase':'6H','inputs':{str(p):sha(p) for p in [P_BIND,P_LEAD,P_AUDIT,P_WALD]},'seed':20261012,'matching':'500 repeated random samples with replacement if needed, same discovery gene_type and global expression quintile; descriptive only','tests':'two-sided Mann-Whitney exploratory, unadjusted, not independent evidence','cautions':['GSE282800 cholesterol-depletion LPPS differs from HAP1 serum-starvation -FBS','bigWig promoter signal is not validated binding or peak calling','GSE282800 binding and discovery GSE267018 share HeLa context, while HAP1 validation is independent','FBS and LPPS signal differences can reflect library normalization and ChIP experimental factors','Overlapping promoters and nearby genes violate independence','Exploratory post-selection comparisons; no mechanistic or causal conclusion']}
(OUT/'RUN_MANIFEST_v1.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('=== PHASE 6H COVERAGE ===');print(coverage.to_string(index=False))
print('=== PHASE 6H OXPHOS GENE SIGNAL ===');print(joined[['gene_symbol','unversioned_id','promoter_eligible','FBS_mean','depletion_mean','depletion_minus_FBS']].to_string(index=False))
print('=== PHASE 6H SIGNAL COMPARISONS (EXPLORATORY) ===');print(pd.DataFrame(rows).to_string(index=False))
print('=== PHASE 6H MATCHED SENSITIVITY ===');print(summary.to_string(index=False))
print('Saved:',OUT)
