#!/usr/bin/env python3
"""Phase 6I: transparent, descriptive prioritization of 13 HAP1 OXPHOS reversal genes.
Run from srebf2-regulatory-signature repository root. Requires pandas and numpy.
No new statistical tests, binding peaks, or causal claims. Existing output is never overwritten.
"""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import pandas as pd

G=Path('results/validation/GSE147770/mitochondrial_robustness_v1/OXPHOS_13_GENES_v1.tsv')
H=Path('results/validation/GSE147770/binding_support_v1/OXPHOS_13_PROMOTER_SIGNAL_v1.tsv')
M=Path('results/validation/GSE147770/binding_support_v1/MATCHED_SENSITIVITY_SUMMARY_v1.tsv')
O=Path('results/validation/GSE147770/candidate_prioritization_v1')
for p in (G,H,M):
    if not p.is_file():sys.exit(f'ERROR missing required file: {p}')
if O.exists():sys.exit(f'ERROR refusing to overwrite existing directory: {O}')
read=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
g,h,m=map(read,(G,H,M))
required_g={'unversioned_id','plus_log2FC','minus_log2FC','interaction_padj','direction'}
required_h={'unversioned_id','FBS_mean','depletion_mean','depletion_minus_FBS','promoter_eligible'}
for name,df,req in [('6G',g,required_g),('6H',h,required_h)]:
    missing=req-set(df.columns)
    if missing:sys.exit(f'ERROR {name} missing columns: {sorted(missing)}')
    if df.unversioned_id.duplicated().any():sys.exit(f'ERROR duplicate {name} gene IDs')
if len(g)!=13 or len(h)!=13 or set(g.unversioned_id)!=set(h.unversioned_id):
    sys.exit('ERROR expected matching sets of 13 unique OXPHOS genes')
# Only transfer binding fields; avoid accidental duplicate columns from Phase 6H's merged table.
fields=['unversioned_id','FBS_mean','depletion_mean','depletion_minus_FBS','promoter_eligible']
d=g.merge(h[fields],on='unversioned_id',validate='one_to_one',how='left')
for c in ['plus_log2FC','minus_log2FC','interaction_padj','FBS_mean','depletion_mean','depletion_minus_FBS']:
    d[c]=pd.to_numeric(d[c],errors='coerce')
if d[list(required_g-{'unversioned_id','direction'})].isna().any().any():sys.exit('ERROR missing key DESeq2 values')
if not d.interaction_padj.between(0,1).all():sys.exit('ERROR invalid interaction adjusted p-values')
if not (d.plus_log2FC.lt(0)&d.minus_log2FC.gt(0)&d.interaction_padj.lt(.05)).all():
    sys.exit('ERROR Phase 6G gene set does not meet expected reversal and FDR criteria')
if d[['FBS_mean','depletion_mean','depletion_minus_FBS']].isna().any().any():sys.exit('ERROR missing promoter signal')
if not np.allclose(d.depletion_mean-d.FBS_mean,d.depletion_minus_FBS,atol=1e-5):
    sys.exit('ERROR inconsistent promoter signal difference')
d['delta_log2FC_minus_vs_plus']=d.minus_log2FC-d.plus_log2FC
d['minus_log10_interaction_padj']=-np.log10(d.interaction_padj.clip(lower=np.finfo(float).tiny))
# Predefined descriptive priority: two equally weighted RNA-only ranks.
# Binding signal is deliberately NOT used for priority due to failed matched-background support.
d['interaction_evidence_percentile']=d.minus_log10_interaction_padj.rank(pct=True,method='average')
d['reversal_magnitude_percentile']=d.delta_log2FC_minus_vs_plus.rank(pct=True,method='average')
d['rna_priority_score_0_100']=50*(d.interaction_evidence_percentile+d.reversal_magnitude_percentile)
d['promoter_signal_nonzero_either']=(d.FBS_mean.gt(0)|d.depletion_mean.gt(0))
d['promoter_signal_increase_depletion']=d.depletion_minus_FBS.gt(0)
d['matched_promoter_enrichment_supported']=False
# These genes were selected for being Hallmark OXPHOS reversal genes; membership is NOT an independent score.
d['hallmark_oxphos_member']=True
d['direct_binding_established']=False
d['functional_oxphos_change_measured']=False
d['evidence_class']='RNA interaction + direction reversal + pathway association; descriptive promoter signal only'
d=d.sort_values(['rna_priority_score_0_100','delta_log2FC_minus_vs_plus','unversioned_id'],ascending=[False,False,True]).reset_index(drop=True)
d.insert(0,'rna_priority_rank',np.arange(1,len(d)+1))
O.mkdir(parents=True)
d.to_csv(O/'CANDIDATE_EVIDENCE_MATRIX_v1.tsv',sep='\t',index=False)
d[['rna_priority_rank','gene_symbol','unversioned_id','rna_priority_score_0_100','interaction_padj','plus_log2FC','minus_log2FC','delta_log2FC_minus_vs_plus','FBS_mean','depletion_mean','depletion_minus_FBS','direct_binding_established']].to_csv(O/'CANDIDATE_SHORTLIST_v1.tsv',sep='\t',index=False)
# Explicitly preserve matched-background caution and audit without pretending gene-level matched tests exist.
m.to_csv(O/'MATCHED_BACKGROUND_CONTEXT_v1.tsv',sep='\t',index=False)
summary=pd.DataFrame([{'genes':len(d),'expected_reversal_genes':int(((d.plus_log2FC<0)&(d.minus_log2FC>0)).sum()),'interaction_FDR05_genes':int((d.interaction_padj<.05).sum()),'nonzero_promoter_signal_genes':int(d.promoter_signal_nonzero_either.sum()),'promoter_increased_under_cholesterol_depletion':int(d.promoter_signal_increase_depletion.sum()),'genes_with_individually_confirmed_binding':0,'functional_oxphos_assays_performed':0}])
summary.to_csv(O/'EVIDENCE_SUMMARY_v1.tsv',sep='\t',index=False)
manifest={'phase':'6I','inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (G,H,M)},'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__,'ranking':'descriptive mean of within-13-gene percentile ranks for -log10(interaction padj) and absolute KO effect difference (minus_FBS - plus_FBS); weights 50/50; NOT a validated prediction model','ranking_scope':'Only 13 post-selected OXPHOS genes, not all candidate genes','promoter_policy':'Promoter signal descriptive only; not ranked because matched-background comparisons did not support consistent enrichment','caveats':['HAP1 -FBS serum starvation differs from HeLa cholesterol depletion','No ChIP peak or binding-specific control validation','Gene-level interaction FDR does not prove significance of each condition-specific fold change','Ranking does not measure direct binding or mitochondrial function','Promoter measurements may reflect technical factors','Pathway membership and gene directions are selection criteria, not independent validation','Matched comparison summary applies to the set, not individual genes']}
(O/'RUN_MANIFEST_v1.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('=== PHASE 6I EVIDENCE SUMMARY ===');print(summary.to_string(index=False))
print('=== PHASE 6I CANDIDATE SHORTLIST (RNA-ONLY PRIORITY) ===')
print(d[['rna_priority_rank','gene_symbol','rna_priority_score_0_100','plus_log2FC','minus_log2FC','delta_log2FC_minus_vs_plus','interaction_padj','FBS_mean','depletion_mean']].to_string(index=False,float_format=lambda x:f'{x:.6g}'))
print('=== PHASE 6I INTERPRETATION ===')
print('Ranking is descriptive and based ONLY on RNA interaction strength and reversal magnitude.')
print('All 13 genes have measured promoter signal, but no direct binding is established.')
print('Matched-background promoter enrichment was not consistently supported in Phase 6H.')
print('Saved:',O)
