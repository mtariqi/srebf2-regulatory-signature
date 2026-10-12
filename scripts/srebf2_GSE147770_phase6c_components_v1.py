#!/usr/bin/env python3
"""Phase 6C: read-only decomposition of frozen Phase 6B HAP1 scores.

Run from the repository root. Requires pandas and numpy. Does not re-score genes.
"""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import pandas as pd

BASE = Path('results/signature_validation/GSE147770_phase6b_v1')
OUT = Path('results/signature_validation/GSE147770_phase6c_v1')
SAMPLE = BASE / 'SAMPLE_SCORES_v1.tsv'
GROUP = BASE / 'GROUP_COMPARISONS_v1.tsv'
METHODS = ['all_gene_percentile_primary','detected_any_percentile','log2_CPM_plus_1','log2_DESeq2_normalized_plus_1']
SERUM = ['plus_FBS','minus_FBS']
GENOTYPES = ['SREBF2_KO','C12orf49_KO']
METRICS = ['positive_component','negative_component','signature_score']

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def main():
    for p in [SAMPLE,GROUP]:
        if not p.is_file(): raise FileNotFoundError(p)
    if OUT.exists(): raise FileExistsError(f'Output directory already exists: {OUT}; do not overwrite.')
    samples=pd.read_csv(SAMPLE,sep='\t')
    groups=pd.read_csv(GROUP,sep='\t')
    required={'sample','genotype','serum','replicate','method',*METRICS}
    if not required.issubset(samples.columns): raise ValueError(f'Missing sample columns: {required-set(samples.columns)}')
    gr_required={'method','serum','comparison','metric','difference'}
    if not gr_required.issubset(groups.columns): raise ValueError(f'Missing comparison columns: {gr_required-set(groups.columns)}')
    if samples[list(METRICS)].isna().any().any(): raise ValueError('Missing score/component values')
    if not np.isfinite(samples[list(METRICS)].to_numpy(dtype=float)).all(): raise ValueError('Nonfinite score/component values')
    if samples.duplicated(['sample','method']).any(): raise ValueError('Duplicate sample/method records')
    if set(samples.method)!=set(METHODS): raise ValueError(f'Unexpected methods: {set(samples.method)^set(METHODS)}')
    if set(samples.serum)!=set(SERUM): raise ValueError('Unexpected serum groups')
    if set(samples.genotype)!={'WT',*GENOTYPES}: raise ValueError('Unexpected genotypes')
    if len(samples)!=72: raise ValueError(f'Expected 72 sample-method records, got {len(samples)}')
    if not np.allclose(samples.positive_component-samples.negative_component,samples.signature_score,rtol=0,atol=1e-10):
        raise ValueError('Signature score != positive minus negative component')
    for method in METHODS:
        for serum in SERUM:
            for genotype in ['WT',*GENOTYPES]:
                x=samples.query('method == @method and serum == @serum and genotype == @genotype')
                if len(x)!=3 or set(x.replicate.astype(int))!={1,2,3}: raise ValueError(f'Unexpected replicates for {method} {serum} {genotype}')
    records=[]; replicate_records=[]; qc=[]
    for method in METHODS:
        for serum in SERUM:
            wt=samples[(samples.method==method)&(samples.serum==serum)&(samples.genotype=='WT')]
            for genotype in GENOTYPES:
                ko=samples[(samples.method==method)&(samples.serum==serum)&(samples.genotype==genotype)]
                d={m:float(ko[m].mean()-wt[m].mean()) for m in METRICS}
                if not np.isclose(d['positive_component']-d['negative_component'],d['signature_score'],atol=1e-10): raise ValueError('Delta decomposition failed')
                comparison=f'{genotype}_minus_WT'
                for metric in METRICS:
                    existing=groups[(groups.method==method)&(groups.serum==serum)&(groups.comparison==comparison)&(groups.metric==metric)]
                    if len(existing)!=1: raise ValueError(f'Cannot verify existing comparison {method} {serum} {comparison} {metric}')
                    diff=float(existing.difference.iloc[0])
                    if not np.isclose(diff,d[metric],atol=1e-9,rtol=1e-9): raise ValueError(f'Phase 6B comparison mismatch: {method} {serum} {comparison} {metric}')
                records.append(dict(method=method,serum=serum,comparison=comparison,n_WT=3,n_perturbed=3,
                    WT_positive_mean=float(wt.positive_component.mean()),perturbed_positive_mean=float(ko.positive_component.mean()),
                    delta_positive=d['positive_component'],WT_negative_mean=float(wt.negative_component.mean()),
                    perturbed_negative_mean=float(ko.negative_component.mean()),delta_negative=d['negative_component'],
                    WT_score_mean=float(wt.signature_score.mean()),perturbed_score_mean=float(ko.signature_score.mean()),
                    delta_score=d['signature_score'],positive_contribution=d['positive_component'],
                    negative_contribution=-d['negative_component'],
                    net_direction='decreased' if d['signature_score']<0 else ('increased' if d['signature_score']>0 else 'unchanged')))
                for _,r in ko.iterrows():
                    rec=dict(method=method,serum=serum,comparison=comparison,sample=r['sample'],replicate=int(r.replicate))
                    for metric in METRICS:
                        rec[metric]=float(r[metric]); rec[f'WT_mean_{metric}']=float(wt[metric].mean())
                        rec[f'delta_vs_WT_mean_{metric}']=float(r[metric]-wt[metric].mean())
                    replicate_records.append(rec)
                qc.append(dict(method=method,serum=serum,comparison=comparison,
                    n_below_WT_mean=int((ko.signature_score<wt.signature_score.mean()).sum()),
                    n_above_WT_mean=int((ko.signature_score>wt.signature_score.mean()).sum()),
                    n_positive_below_WT_mean=int((ko.positive_component<wt.positive_component.mean()).sum()),
                    n_negative_below_WT_mean=int((ko.negative_component<wt.negative_component.mean()).sum())))
    out=pd.DataFrame(records); rep=pd.DataFrame(replicate_records); consistency=pd.DataFrame(qc)
    interaction=[]
    for method in METHODS:
        for genotype in GENOTYPES:
            sub=out[(out.method==method)&(out.comparison==f'{genotype}_minus_WT')].set_index('serum')
            row={'method':method,'comparison':f'{genotype}_minus_WT'}
            for metric,col in [('positive','delta_positive'),('negative','delta_negative'),('score','delta_score')]:
                row[f'{metric}_plus_FBS']=float(sub.loc['plus_FBS',col]); row[f'{metric}_minus_FBS']=float(sub.loc['minus_FBS',col]);
                row[f'{metric}_difference_of_differences']=row[f'{metric}_minus_FBS']-row[f'{metric}_plus_FBS']
            interaction.append(row)
    interactions=pd.DataFrame(interaction)
    OUT.mkdir(parents=True)
    out.to_csv(OUT/'COMPONENT_DECOMPOSITION_v1.tsv',sep='\t',index=False,float_format='%.12g')
    rep.to_csv(OUT/'REPLICATE_COMPONENTS_v1.tsv',sep='\t',index=False,float_format='%.12g')
    consistency.to_csv(OUT/'REPLICATE_CONSISTENCY_v1.tsv',sep='\t',index=False)
    interactions.to_csv(OUT/'SERUM_DIFFERENCE_OF_DIFFERENCES_v1.tsv',sep='\t',index=False,float_format='%.12g')
    manifest={'phase':'6C','version':'v1','analysis':'descriptive decomposition of frozen Phase 6B scores',
        'source_files_sha256':{str(SAMPLE):sha256(SAMPLE),str(GROUP):sha256(GROUP)},
        'methods':METHODS,'serum_conditions':SERUM,'perturbations':GENOTYPES,
        'definition':'signature = positive - negative; component deltas = perturbed mean - WT mean; serum difference of differences = minus_FBS delta - plus_FBS delta',
        'no_rescoring':True,'no_inferential_p_values':True,'replicates_per_group':3}
    (OUT/'RUN_MANIFEST_v1.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('=== SREBF2 KO COMPONENT DECOMPOSITION ===')
    print(out[out.comparison=='SREBF2_KO_minus_WT'][['method','serum','delta_positive','delta_negative','delta_score']].to_string(index=False,float_format=lambda v:f'{v:+.6f}'))
    print('\n=== SREBF2 KO REPLICATE CONSISTENCY ===')
    print(consistency[consistency.comparison=='SREBF2_KO_minus_WT'].to_string(index=False))
    print('\n=== SREBF2 KO SERUM DIFFERENCE OF DIFFERENCES (DESCRIPTIVE) ===')
    print(interactions[interactions.comparison=='SREBF2_KO_minus_WT'][['method','positive_difference_of_differences','negative_difference_of_differences','score_difference_of_differences']].to_string(index=False,float_format=lambda v:f'{v:+.6f}'))
    print('\nSaved:',OUT.resolve())

if __name__=='__main__':
    try: main()
    except Exception as exc:
        print(f'ERROR: {exc}',file=sys.stderr); sys.exit(1)
