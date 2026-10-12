#!/usr/bin/env python3
"""GSE271000 Phase 4C: common-support background matching (descriptive)."""
from pathlib import Path
import gzip, hashlib, json
import numpy as np
import pandas as pd

ROOT=Path.cwd()
AUDIT=ROOT/'results/signature/discovery_v1/SIGNATURE_SELECTION_AUDIT.tsv'
SIG=ROOT/'results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv'
RAW=ROOT/'data/raw/validation/GSE271000'
OUT=ROOT/'results/validation/GSE271000/common_support_v3'
SEED=20261011; DRAWS=200
if OUT.exists(): raise SystemExit(f'Refusing overwrite: {OUT}')
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def yes(x):return x.astype(str).str.lower().eq('true')
a=pd.read_csv(AUDIT,sep='\t',dtype={'gene_id':str,'gene_symbol':str})
s=pd.read_csv(SIG,sep='\t',dtype={'gene_id':str,'gene_symbol':str})
assert len(a)==17158 and len(s)==1326
assert a.gene_id.is_unique and s.gene_id.is_unique
assert set(a.loc[yes(a.signature_selected),'gene_id'])==set(s.gene_id)
a['selected']=yes(a.signature_selected)
a['lfc']=pd.to_numeric(a.log2FC_primary,errors='coerce')
a['base']=pd.to_numeric(a.baseMean,errors='coerce')
a['biotype']=a.gene_type.fillna('unknown').astype(str)
a['biotype_broad']=np.where(a.biotype.eq('protein_coding'),'protein_coding','other')
a['log_expr']=np.log1p(a.base.clip(lower=0))
a['symbol_unique']=a.gene_symbol.notna() & a.gene_symbol.ne('') & a.gene_symbol.map(a.gene_symbol.value_counts()).eq(1)
a['expected_sign']=np.sign(a.lfc)
a['direction_group']=np.where(a.selected,a.direction,'unselected')
# Define all strata from discovery data, independent of validation outcomes.
a['decile']=pd.qcut(a.log_expr.rank(method='first'),10,labels=False)
a['quintile']=pd.qcut(a.log_expr.rank(method='first'),5,labels=False)
a['strict']=a.biotype+'|'+a.decile.astype(str)
a['broad']=a.biotype_broad+'|'+a.quintile.astype(str)
rng=np.random.default_rng(SEED)
feas=[]; results=[]; draws=[]; coverage=[]; hashes={str(AUDIT):digest(AUDIT),str(SIG):digest(SIG)}
OUT.mkdir(parents=True)
for hp in ('664','665','667'):
    p=RAW/f'GSE271000_DESeq_HepG2_shSREBP2_{hp}.txt.gz';hashes[str(p)]=digest(p)
    with gzip.open(p,'rt') as f:d=pd.read_csv(f,sep='\t',index_col=0)
    d.index=d.index.astype(str)
    duplicates=set(d.index[d.index.duplicated(keep=False)])
    d=d.loc[~d.index.isin(duplicates)]
    assert 'log2FoldChange' in d
    z=a.copy();z['validation_lfc']=z.gene_symbol.map(pd.to_numeric(d.log2FoldChange,errors='coerce'))
    z.loc[~z.symbol_unique | z.gene_symbol.isin(duplicates),'validation_lfc']=np.nan
    z['evaluable']=z.validation_lfc.notna() & np.isfinite(z.validation_lfc) & z.validation_lfc.ne(0) & z.expected_sign.ne(0)
    z['agree']=z.expected_sign.eq(np.sign(z.validation_lfc))
    e=z.loc[z.evaluable]
    coverage.append({'hairpin':hp,'signature_total':int(z.selected.sum()),'signature_evaluable':int((z.selected&z.evaluable).sum()),'background_evaluable':int((~z.selected&z.evaluable).sum())})
    bg=e.loc[~e.selected]
    for group in ('decreased_in_KO','increased_in_KO','all_selected'):
        case=e.loc[e.selected & ((e.direction_group==group) if group!='all_selected' else True)]
        for scheme,col in [('strict_biotype_decile','strict'),('broad_biotype_quintile','broad')]:
            cc=case[col].value_counts();bc=bg[col].value_counts()
            for st,n in cc.items():
                available=int(bc.get(st,0))
                feas.append({'hairpin':hp,'group':group,'scheme':scheme,'stratum':st,'case_n':int(n),'control_n':available,'shortfall':max(0,int(n)-available)})
            # Restrict to supported strata, and within each retain at most N_controls
            # cases, chosen without using validation concordance outcomes.
            kept=[]; pools={}; dropped=0
            for st,n in cc.items():
                eligible=case.index[case[col]==st].to_numpy()
                pool=bg.index[bg[col]==st].to_numpy()
                take=min(len(eligible),len(pool))
                if take:
                    # Deterministic stable order by gene ID, not validation outcome.
                    chosen=case.loc[eligible].sort_values('gene_id').index.to_numpy()[:take]
                    kept.extend(chosen)
                    pools[st]=(pool,take)
                dropped+=len(eligible)-take
            sub=case.loc[kept]
            rec={'hairpin':hp,'group':group,'scheme':scheme,'total_evaluable_case':len(case),
                 'common_support_case':len(sub),'excluded_case':dropped,
                 'retention_fraction':len(sub)/len(case) if len(case) else np.nan,
                 'case_concordance_full':case.agree.mean(),
                 'case_concordance_supported':sub.agree.mean() if len(sub) else np.nan,
                 'matching_status':'matched' if len(sub) else 'no_support'}
            results.append(rec)
            if not len(sub):continue
            for draw in range(1,DRAWS+1):
                chosen=np.concatenate([rng.choice(pool,size=k,replace=False) for pool,k in pools.values()])
                control=bg.loc[chosen]
                draws.append({'hairpin':hp,'group':group,'scheme':scheme,'draw':draw,
                              'case_n':len(sub),'control_n':len(control),
                              'case_concordance':float(sub.agree.mean()),
                              'control_concordance':float(control.agree.mean()),
                              'difference':float(sub.agree.mean()-control.agree.mean())})
    z[['gene_id','gene_symbol','selected','direction_group','biotype','base','lfc','strict','broad','evaluable','validation_lfc','agree']].to_csv(OUT/f'GENE_AUDIT_{hp}_v3.tsv',sep='\t',index=False)
feas=pd.DataFrame(feas);res=pd.DataFrame(results);dr=pd.DataFrame(draws)
feas.to_csv(OUT/'STRATUM_FEASIBILITY_v3.tsv',sep='\t',index=False)
res.to_csv(OUT/'COMMON_SUPPORT_SUMMARY_v3.tsv',sep='\t',index=False)
dr.to_csv(OUT/'MATCHED_DRAWS_v3.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(OUT/'COVERAGE_v3.tsv',sep='\t',index=False)
manifest={'version':'common_support_v3','seed':SEED,'draws':DRAWS,'sha256':hashes,
          'strata':{'strict':'exact biotype x global discovery expression decile','broad':'protein_coding/other x global discovery expression quintile'},
          'case_selection':'For each stratum retain first min(case_n,control_n) cases sorted by gene_id; independent of validation LFC; excluded cases reported.',
          'draws':'Control sampling without replacement within each draw; draws not independent CIs.',
          'limitations':['Common-support results only apply to retained genes; no extrapolation to all 1326.',
                         'Gene sign correlations and small discovery LFC among controls make this descriptive, not a calibrated null.',
                         'Assumes shSREBP2 versus shControl orientation; independently verify metadata.',
                         'No genome-block inference or independent gene assumption.']}
(OUT/'COMMON_SUPPORT_MANIFEST_v3.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('\nCoverage:\n',pd.DataFrame(coverage).to_string(index=False))
print('\nCommon support:\n',res.to_string(index=False))
print('\nMatched differences across 200 draws (not CIs):')
if len(dr):print(dr.groupby(['hairpin','group','scheme']).difference.agg(['min','median','max']).to_string())
print('\nSaved:',OUT)
