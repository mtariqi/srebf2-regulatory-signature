#!/usr/bin/env python3
"""Phase 6B external HAP1 scoring. Frozen gene directions, cohort-adapted mapping/universe.
Requires pandas/numpy. Never overwrites outputs. Run from repository root.
"""
import hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path.cwd()
COUNT=ROOT/'data/raw/validation/GSE147770/GSE147770_STAR_readCounts.txt.gz'
SIG=ROOT/'results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv'
AUDIT=ROOT/'results/validation/GSE147770/audit_v4'
OUT=ROOT/'results/signature_validation/GSE147770_phase6b_v1'
SAMPLES=[f'{g}{i}{s}' for g in 'WSC' for i in range(1,4) for s in ('plus','min')]
EXPECTED='3f4a2307ede59793b77ac50c38fa519ece012c0a2128fe56ea4c0f44f0b7d598'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def main():
    for p in [COUNT,SIG,AUDIT/'MATRIX_ROW_AUDIT_v4.tsv',AUDIT/'SIGNATURE_MAPPING_v4.tsv',AUDIT/'AUDIT_MANIFEST_v4.json']:
        if not p.is_file():raise FileNotFoundError(p)
    if OUT.exists():raise FileExistsError(f'Refusing to overwrite: {OUT}')
    if sha(COUNT)!=EXPECTED:raise ValueError('Count SHA256 mismatch')
    manifest=json.loads((AUDIT/'AUDIT_MANIFEST_v4.json').read_text())
    if manifest['count_sha256']!=sha(COUNT) or manifest['signature_sha256']!=sha(SIG):raise ValueError('Audit/source checksum mismatch')
    m=pd.read_csv(COUNT,sep='\t',dtype={'Ensembl.ID':'string','Gene.symbol':'string','Gene.ID':'string'},low_memory=False)
    audit=pd.read_csv(AUDIT/'MATRIX_ROW_AUDIT_v4.tsv',sep='\t')
    mapping=pd.read_csv(AUDIT/'SIGNATURE_MAPPING_v4.tsv',sep='\t',dtype='string')
    sig=pd.read_csv(SIG,sep='\t',dtype={'gene_id':'string'})
    if len(m)!=38344 or len(audit)!=len(m) or len(sig)!=1326:raise ValueError('Unexpected row count')
    if set(m.columns.tolist()[3:])!=set(SAMPLES) or len(m.columns.tolist()[3:])!=18:raise ValueError('Unexpected sample columns')
    if not np.array_equal(audit['Ensembl.ID'].fillna('').astype(str),m['Ensembl.ID'].fillna('').astype(str)):
        raise ValueError('Audit matrix row mismatch')
    if not mapping['gene_id'].astype(str).equals(sig['gene_id'].astype(str)):
        raise ValueError('Mapping/signature row mismatch')
    eligible=audit['eligible_complete_unique_gene'].astype(bool).to_numpy()
    if eligible.sum()!=37919:raise ValueError('Unexpected eligible gene universe')
    raw=m.loc[eligible,SAMPLES].apply(pd.to_numeric,errors='raise').to_numpy(dtype=np.float64)
    if not np.isfinite(raw).all() or (raw<0).any() or not np.equal(raw,np.floor(raw)).all():
        raise ValueError('Eligible matrix has invalid counts')
    ids=m.loc[eligible,'Ensembl.ID'].astype(str).str.replace(r'\.\d+$','',regex=True).tolist()
    if len(set(ids))!=len(ids):raise ValueError('Eligible IDs duplicated')
    rowmap={v:i for i,v in enumerate(ids)}
    sel=mapping['mapping_status'].eq('unique_complete').to_numpy()
    if sel.sum()!=1065:raise ValueError('Unexpected signature coverage')
    gene_ids=mapping.loc[sel,'unversioned_id'].astype(str).tolist()
    ix=np.array([rowmap[g] for g in gene_ids],dtype=int)
    directions=mapping.loc[sel,'direction'].astype(str).to_numpy()
    pos=ix[directions=='decreased_in_KO'];neg=ix[directions=='increased_in_KO']
    if len(pos)!=676 or len(neg)!=389:raise ValueError('Unexpected directional mapping')
    # DESeq2-style size factors: only genes strictly positive across all 18 samples.
    shared=(raw>0).all(axis=1)
    if not shared.any():raise ValueError('No all-positive genes for normalization')
    loggeom=np.log(raw[shared]).mean(axis=1)
    size=np.median(np.exp(np.log(raw[shared])-loggeom[:,None]),axis=0)
    if not np.isfinite(size).all() or (size<=0).any():raise ValueError('Invalid size factors')
    norm=raw/size[None,:]
    # Ranking midranks, 1-based: (rank - 0.5)/N, including zero-count genes.
    pct=pd.DataFrame(norm).rank(axis=0,method='average',pct=False).to_numpy()
    pct=(pct-0.5)/len(raw)
    detected=(raw>0).any(axis=1)
    if detected.sum()==0:raise ValueError('No detected genes')
    pct_detected=pd.DataFrame(norm[detected]).rank(axis=0,method='average').to_numpy()
    pct_detected=(pct_detected-0.5)/detected.sum()
    detected_indices=np.flatnonzero(detected)
    inv=np.full(len(raw),-1,dtype=int);inv[detected_indices]=np.arange(len(detected_indices))
    # Signature genes with all-zero counts have no detected-universe rank; omit them only for this sensitivity.
    logcpm=np.log2(raw/raw.sum(axis=0)[None,:]*1e6+1)
    lognorm=np.log2(norm+1)
    methods={
        'all_gene_percentile_primary':(pct,pos,neg),
        'detected_any_percentile':(pct_detected,inv[pos][inv[pos]>=0],inv[neg][inv[neg]>=0]),
        'log2_CPM_plus_1':(logcpm,pos,neg),
        'log2_DESeq2_normalized_plus_1':(lognorm,pos,neg),
    }
    records=[]
    for method,(mat,pi,ni) in methods.items():
        if len(pi)==0 or len(ni)==0:raise ValueError(f'Empty signature component: {method}')
        pv=mat[pi].mean(axis=0);nv=mat[ni].mean(axis=0)
        for j,s in enumerate(SAMPLES):
            records.append(dict(sample=s,genotype={'W':'WT','S':'SREBF2_KO','C':'C12orf49_KO'}[s[0]],serum='plus_FBS' if s.endswith('plus') else 'minus_FBS',replicate=int(s[1]),method=method,positive_component=float(pv[j]),negative_component=float(nv[j]),signature_score=float(pv[j]-nv[j]),size_factor=float(size[j]),eligible_library_size=int(raw[:,j].sum()),detected_genes=int((raw[:,j]>0).sum()),positive_genes_scored=len(pi),negative_genes_scored=len(ni)))
    scores=pd.DataFrame(records)
    summaries=[]
    for (method,serum),sub in scores.groupby(['method','serum'],sort=False):
        wt=sub[sub.genotype=='WT']
        for genotype in ['SREBF2_KO','C12orf49_KO']:
            ko=sub[sub.genotype==genotype]
            if len(wt)!=3 or len(ko)!=3:raise ValueError('Group replicate mismatch')
            for metric in ['signature_score','positive_component','negative_component']:
                summaries.append(dict(method=method,serum=serum,comparison=f'{genotype}_minus_WT',metric=metric,n_WT=len(wt),n_perturbed=len(ko),WT_mean=float(wt[metric].mean()),perturbed_mean=float(ko[metric].mean()),difference=float(ko[metric].mean()-wt[metric].mean()),all_perturbed_below_WT_mean=bool((ko[metric]<wt[metric].mean()).all())))
    comparisons=pd.DataFrame(summaries)
    # SREBF2 index expression is an independent, descriptive benchmark, not part of the signature.
    index_rows=[i for i,g in enumerate(ids) if g=='ENSG00000072310']  # SREBF2 Ensembl stable ID
    index_data=[]
    if len(index_rows)==1:
        i=index_rows[0]
        for j,s in enumerate(SAMPLES):
            index_data.append(dict(sample=s,genotype={'W':'WT','S':'SREBF2_KO','C':'C12orf49_KO'}[s[0]],serum='plus_FBS' if s.endswith('plus') else 'minus_FBS',SREBF2_raw_count=int(raw[i,j]),SREBF2_log2_normalized_plus_1=float(lognorm[i,j])))
    else:
        print('WARNING: SREBF2 stable ID not uniquely found; index benchmark unavailable')
    OUT.mkdir(parents=True,exist_ok=False)
    scores.to_csv(OUT/'SAMPLE_SCORES_v1.tsv',sep='\t',index=False)
    comparisons.to_csv(OUT/'GROUP_COMPARISONS_v1.tsv',sep='\t',index=False)
    pd.DataFrame(index_data,columns=['sample','genotype','serum','SREBF2_raw_count','SREBF2_log2_normalized_plus_1']).to_csv(OUT/'SREBF2_INDEX_EXPRESSION_v1.tsv',sep='\t',index=False)
    output_manifest={'phase':'6B','status':'external_cohort_adaptation_of_frozen_phase5_policy','source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [COUNT,SIG,AUDIT/'MATRIX_ROW_AUDIT_v4.tsv',AUDIT/'SIGNATURE_MAPPING_v4.tsv',ROOT/'metadata/SREBF2_VALIDATION_SCORING_POLICY_v1.json']},'eligible_genes':int(len(raw)),'signature_positive_mapped':int(len(pos)),'signature_negative_mapped':int(len(neg)),'normalization_all_positive_genes':int(shared.sum()),'detected_any_genes':int(detected.sum()),'index_gene_found':len(index_rows)==1,'methods':list(methods),'no_source_modification':True,'notes':['HAP1-specific eligible universe and gene mapping; not a literal application of Phase 5 exact-ID policy','No inferential significance testing; comparisons are descriptive','C12orf49 is exploratory, not WT','Different datasets scores are not directly comparable']}
    (OUT/'RUN_MANIFEST_v1.json').write_text(json.dumps(output_manifest,indent=2)+'\n')
    print('=== PHASE 6B SUMMARY ===')
    print('Eligible genes:',len(raw),'Signature positive:',len(pos),'negative:',len(neg),'Normalization genes:',int(shared.sum()))
    print('=== SREBF2 KO MINUS WT: SIGNATURE SCORES ===')
    print(comparisons[(comparisons.comparison=='SREBF2_KO_minus_WT')&(comparisons.metric=='signature_score')][['method','serum','WT_mean','perturbed_mean','difference','all_perturbed_below_WT_mean']].to_string(index=False))
    print('=== SREBF2 INDEX EXPRESSION ===')
    if index_data:print(pd.DataFrame(index_data).to_string(index=False))
    print('Saved:',OUT)

if __name__=='__main__':main()
