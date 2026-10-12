#!/usr/bin/env python3
"""Phase 5D exploratory benchmark and expression/biotype-matched gene-set controls.

Run from repository root. Read-only inputs; refuses to overwrite outputs.
Requires pandas/numpy. No inferential p-values or cross-validation claims.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

ROOT=Path.cwd()
COUNT=ROOT/'results/counts/validation/GSE271000/GSE271000_validation_gene_counts.tsv'
META=ROOT/'metadata/GSE271000_validation_samples.tsv'
POLICY=ROOT/'metadata/SREBF2_VALIDATION_SCORING_POLICY_v1.json'
SIG=ROOT/'results/signature/discovery_v1'
AUDIT=SIG/'SIGNATURE_SELECTION_AUDIT.tsv'
SAMPLES=ROOT/'results/signature_validation/GSE271000_sensitivity_v1/SAMPLE_METHOD_SCORES_v1.tsv'
SUMMARY=ROOT/'results/signature_validation/GSE271000_sensitivity_v1/GROUP_METHOD_COMPARISON_v1.tsv'
OUT=ROOT/'results/signature_validation/GSE271000_benchmark_v2'
SEED=20261011
DRAWS=200
SREBF2='ENSG00000198911.13'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def main():
    if OUT.exists():raise SystemExit(f'Refusing to overwrite existing output: {OUT}')
    required=[COUNT,META,POLICY,AUDIT,SAMPLES,SUMMARY,SIG/'SIGNATURE_MANIFEST.json',SIG/'SREBF2_SIGNATURE_POSITIVE.tsv',SIG/'SREBF2_SIGNATURE_NEGATIVE.tsv']
    for f in required:
        if not f.is_file():raise SystemExit(f'Missing required file: {f}')
    policy=json.loads(POLICY.read_text())
    for f,k in [(SIG/'SIGNATURE_MANIFEST.json','discovery_manifest_sha256'),(SIG/'SREBF2_SIGNATURE_POSITIVE.tsv','positive_signature_sha256'),(SIG/'SREBF2_SIGNATURE_NEGATIVE.tsv','negative_signature_sha256')]:
        if sha(f)!=policy[k]:raise ValueError(f'Frozen policy signature checksum mismatch: {f}')
    meta=pd.read_csv(META,sep='\t',dtype=str)
    expected={'control':2,'shSREBP2_664':2,'shSREBP2_665':2,'shSREBP2_667':2}
    if len(meta)!=8 or meta.run.duplicated().any() or meta.groupby('group').size().to_dict()!=expected:raise ValueError('Unexpected sample metadata')
    data=pd.read_csv(COUNT,sep='\t',dtype={'gene_id':str})
    if data.shape!=(78932,9) or data.columns.tolist()!=['gene_id']+meta.run.tolist() or data.gene_id.duplicated().any():raise ValueError('Unexpected count matrix')
    counts=data.iloc[:,1:].to_numpy(dtype=np.float64)
    if not np.isfinite(counts).all() or np.any(counts<0) or np.any(counts!=np.floor(counts)):raise ValueError('Invalid count values')
    genes=pd.Index(data.gene_id)
    pos=pd.read_csv(SIG/'SREBF2_SIGNATURE_POSITIVE.tsv',sep='\t',dtype={'gene_id':str}).gene_id.tolist()
    neg=pd.read_csv(SIG/'SREBF2_SIGNATURE_NEGATIVE.tsv',sep='\t',dtype={'gene_id':str}).gene_id.tolist()
    if len(pos)!=795 or len(neg)!=531 or len(set(pos+neg))!=1326:raise ValueError('Frozen signature size mismatch')
    if any(x not in genes for x in pos+neg+[SREBF2]):raise ValueError('Missing signature or SREBF2 gene in counts')
    audit=pd.read_csv(AUDIT,sep='\t',dtype={'gene_id':str})
    if len(audit)!=17158 or audit.gene_id.duplicated().any():raise ValueError('Unexpected discovery audit')
    if not {'gene_id','gene_type','baseMean','signature_selected'}.issubset(audit.columns):raise ValueError('Missing audit columns')
    audit=audit[audit.gene_id.isin(genes)].copy()
    audit['selected']=audit.gene_id.isin(set(pos+neg))
    if audit.selected.sum()!=1326:raise ValueError('Audit/signature selection mismatch')
    audit['biotype2']=np.where(audit.gene_type.eq('protein_coding'),'protein_coding','other')
    audit['expr']=np.log1p(pd.to_numeric(audit.baseMean,errors='coerce').clip(lower=0))
    if audit.expr.isna().any():raise ValueError('Missing discovery expression for matched pool')
    audit['expr_bin']=pd.qcut(audit.expr.rank(method='first'),q=5,labels=False,duplicates='drop')
    audit['stratum']=audit.biotype2+'|'+audit.expr_bin.astype(str)
    background=audit[~audit.selected & ~audit.gene_id.eq(SREBF2)]
    sample=pd.read_csv(SAMPLES,sep='\t')
    if set(sample.run)!=set(meta.run) or set(sample.group)!=set(meta.group):raise ValueError('Sensitivity sample mismatch')
    methods=['primary_all_gene_percentile','detected_any_gene_percentile','log2_CPM_plus1','log2_DESeq2_normalized_plus1']
    if len(sample)!=32 or set(sample.method)!=set(methods):raise ValueError('Unexpected sensitivity method table')
    # Recompute the same normalization arrays as Phase 5C for gene-set controls.
    valid=(counts>0).all(axis=1)
    geom=np.exp(np.log(counts[valid]).mean(axis=1))
    sf=np.median(counts[valid]/geom[:,None],axis=0)
    library=counts.sum(axis=0)
    original=(pd.DataFrame(counts).rank(axis=0,method='average').to_numpy()-.5)/len(genes)
    detected=(counts>0).any(axis=1)
    detected_mat=np.full(counts.shape,np.nan)
    detected_mat[detected]=(pd.DataFrame(counts[detected]).rank(axis=0,method='average').to_numpy()-.5)/int(detected.sum())
    matrices={'primary_all_gene_percentile':original,'detected_any_gene_percentile':detected_mat,'log2_CPM_plus1':np.log2(counts/library[None,:]*1e6+1),'log2_DESeq2_normalized_plus1':np.log2(counts/sf[None,:]+1)}
    # Ensure null method matrices reproduce all frozen sensitivity scores.
    pi=genes.get_indexer(pos);ni=genes.get_indexer(neg);si=genes.get_loc(SREBF2)
    for method,mat in matrices.items():
        sub=sample[sample.method.eq(method)].set_index('run').loc[meta.run]
        observed=np.nanmean(mat[pi],axis=0)-np.nanmean(mat[ni],axis=0)
        if not np.allclose(observed,sub.signature_score.to_numpy(float),atol=1e-8):raise ValueError(f'Cannot reproduce {method} score')
        if not np.allclose(mat[si],sub.SREBF2_value.to_numpy(float),atol=1e-8):raise ValueError(f'Cannot reproduce {method} SREBF2')
    # Descriptive metrics on the same eight samples. Standardization uses pooled within-group SD;
    # no p-values, predictive AUC, or cross-validation with n=2 per group.
    group_ix={g:np.flatnonzero(meta.group.eq(g).to_numpy()) for g in expected}
    ctl=group_ix['control']
    metrics=[]
    for method in methods:
        sub=sample[sample.method.eq(method)].set_index('run').loc[meta.run]
        for metric,col in [('full_signature','signature_score'),('positive_component','positive_component'),('negative_component','negative_component'),('SREBF2_expression','SREBF2_value')]:
            values=sub[col].to_numpy(float)
            for hp in ['664','665','667']:
                ix=group_ix['shSREBP2_'+hp]
                d=values[ix].mean()-values[ctl].mean()
                sd=np.sqrt((values[ix].var(ddof=1)+values[ctl].var(ddof=1))/2)
                metrics.append({'method':method,'metric':metric,'hairpin':hp,'delta_vs_control':d,'within_group_pooled_sd':sd,'standardized_delta':d/sd if sd>0 else np.nan,'both_replicates_below_control_mean':bool(np.all(values[ix]<values[ctl].mean())),'both_replicates_above_control_mean':bool(np.all(values[ix]>values[ctl].mean()))})
    metrics=pd.DataFrame(metrics)
    # Gene-set null uses discovery biotype and expression quintile matching, no reuse within each draw.
    # The null is descriptive: it does not account for selection on discovery effect size or gene correlation.
    audit_by_id=audit.set_index('gene_id')
    pool={st:background.loc[background.stratum.eq(st),'gene_id'].to_numpy() for st in background.stratum.unique()}
    demand=pd.concat([audit_by_id.loc[pos,'stratum'],audit_by_id.loc[neg,'stratum']]).value_counts()
    feasible=[]
    for st,n in demand.items():
        available=len(pool.get(st,[]))
        feasible.append({'stratum':st,'signature_n':int(n),'background_n':available,'sufficient':available>=n})
    feasible=pd.DataFrame(feasible)
    if not feasible.sufficient.all():
        # Preserve audit, but do not silently use replacement or alter the signature.
        OUT.mkdir(parents=True,exist_ok=False)
        feasible.to_csv(OUT/'MATCHING_FEASIBILITY_v2.tsv',sep='\t',index=False)
        (OUT/'STOP_REASON.txt').write_text('Matched null could not be run without replacement. Inspect MATCHING_FEASIBILITY_v2.tsv. No benchmark/null results generated.\n')
        raise SystemExit('Insufficient background for exact biotype2 x discovery-expression quintile matching; see output feasibility audit. No overwrite or replacement sampling performed.')
    rng=np.random.default_rng(SEED)
    pstr=audit_by_id.loc[pos,'stratum'].value_counts().to_dict()
    nstr=audit_by_id.loc[neg,'stratum'].value_counts().to_dict()
    draws=[]
    for draw in range(1,DRAWS+1):
        chosen_pos=[];chosen_neg=[]
        for st in sorted(demand.index):
            ids=rng.choice(pool[st],size=int(demand[st]),replace=False)
            n_p=pstr.get(st,0)
            chosen_pos.extend(ids[:n_p]);chosen_neg.extend(ids[n_p:])
        idxp=genes.get_indexer(chosen_pos);idxn=genes.get_indexer(chosen_neg)
        for method,mat in matrices.items():
            score=np.nanmean(mat[idxp],axis=0)-np.nanmean(mat[idxn],axis=0)
            for hp in ['664','665','667']:
                ix=group_ix['shSREBP2_'+hp]
                draws.append({'draw':draw,'method':method,'hairpin':hp,'null_delta':float(score[ix].mean()-score[ctl].mean()),'null_both_replicates_below':bool(np.all(score[ix]<score[ctl].mean()))})
    draws=pd.DataFrame(draws)
    comp=[]
    for (method,hp),g in draws.groupby(['method','hairpin'],sort=False):
        obs=float(metrics[(metrics.method==method)&(metrics.metric=='full_signature')&(metrics.hairpin==hp)].delta_vs_control.iloc[0])
        # Signed negative deltas are expected. This is a descriptive empirical tail fraction,
        # NOT a calibrated significance test because matching does not address all confounding.
        comp.append({'method':method,'hairpin':hp,'observed_signature_delta':obs,'null_delta_median':g.null_delta.median(),'null_delta_min':g.null_delta.min(),'null_delta_max':g.null_delta.max(),'null_fraction_at_least_as_negative_as_observed':float((g.null_delta<=obs).mean()),'null_fraction_both_replicates_below':g.null_both_replicates_below.mean()})
    comp=pd.DataFrame(comp)
    OUT.mkdir(parents=True,exist_ok=False)
    metrics.to_csv(OUT/'BENCHMARK_METRICS_v2.tsv',sep='\t',index=False,float_format='%.10g')
    feasible.to_csv(OUT/'MATCHING_FEASIBILITY_v2.tsv',sep='\t',index=False)
    draws.to_csv(OUT/'MATCHED_NULL_DRAWS_v2.tsv',sep='\t',index=False,float_format='%.10g')
    comp.to_csv(OUT/'MATCHED_NULL_COMPARISON_v2.tsv',sep='\t',index=False,float_format='%.10g')
    manifest={'phase':'5D exploratory benchmarking v2 (NaN-aware reproduction fix)','seed':SEED,'draws':DRAWS,'matching':'discovery gene_type collapsed to protein_coding/other x expression quintile, without replacement; preserve 795/531 sign group sizes','input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in required},'limitations':['Eight samples only: two per group; standardized effects unstable','Three shRNAs share a common control and are not independent cohorts','Alternative scoring methods were explored after viewing validation results','Matched random gene sets do not account for discovery effect-size selection, pathway correlation or latent biology','Random-draw empirical fractions are descriptive, not calibrated p-values','Detected-gene percentile means omit undetected-universe NaNs exactly as Phase 5C; matched gene-set denominators may vary by detected coverage','No claim of binding causality, classification generalization, or superiority over SREBF2 expression']}
    (OUT/'BENCHMARK_MANIFEST_v2.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('=== STANDARDIZED SEPARATION (descriptive; same eight samples) ===')
    print(metrics[metrics.metric.isin(['full_signature','SREBF2_expression'])][['method','metric','hairpin','delta_vs_control','standardized_delta','both_replicates_below_control_mean']].to_string(index=False,float_format=lambda x:f'{x:+.5f}'))
    print('\n=== MATCHED RANDOM GENE-SET CONTROL (descriptive, not p-values) ===')
    print(comp.to_string(index=False,float_format=lambda x:f'{x:+.5f}'))
    print('\nSaved:',OUT)

if __name__=='__main__':main()
