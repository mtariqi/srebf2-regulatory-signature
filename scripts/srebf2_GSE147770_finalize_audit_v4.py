#!/usr/bin/env python3
"""GSE147770 Phase 6A: non-destructive complete-case and signature mapping audit."""
import hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path.cwd()
COUNTS=ROOT/'data/raw/validation/GSE147770/GSE147770_STAR_readCounts.txt.gz'
SIG=ROOT/'results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv'
OUT=ROOT/'results/validation/GSE147770/audit_v4'
SAMPLES=[f'{g}{i}{s}' for g in 'WSC' for i in range(1,4) for s in ('plus','min')]
EXPECTED='3f4a2307ede59793b77ac50c38fa519ece012c0a2128fe56ea4c0f44f0b7d598'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 for p in (COUNTS,SIG):
  if not p.is_file():raise FileNotFoundError(p)
 if OUT.exists() and any(OUT.iterdir()):raise FileExistsError(f'Output exists: {OUT}')
 digest=sha(COUNTS)
 if digest!=EXPECTED:raise ValueError('Count matrix checksum mismatch')
 m=pd.read_csv(COUNTS,sep='\t',dtype={'Ensembl.ID':'string','Gene.symbol':'string','Gene.ID':'string'},low_memory=False)
 if len(m)!=38344:raise ValueError(f'Unexpected row count {len(m)}')
 if any(c not in m.columns for c in SAMPLES):raise ValueError('Missing sample columns')
 if len([c for c in m.columns if re.fullmatch(r'[WSC][123](?:plus|min)',c)])!=18:raise ValueError('Unexpected sample layout')
 vals=m[SAMPLES].apply(pd.to_numeric,errors='coerce')
 invalid=~np.isfinite(vals.to_numpy(dtype=float)) | (vals.to_numpy(dtype=float)<0) | (vals.to_numpy(dtype=float)%1!=0)
 anybad=invalid.any(axis=1); allbad=invalid.all(axis=1); partial=anybad & ~allbad
 ids=m['Ensembl.ID'].astype('string').str.replace(r'\.\d+$','',regex=True)
 standard=ids.str.fullmatch(r'ENSG\d{11}').fillna(False)
 duplicate=ids.duplicated(keep=False) & ids.notna()
 eligible=(~anybad)&standard.to_numpy(dtype=bool)&(~duplicate.to_numpy(dtype=bool))
 sig=pd.read_csv(SIG,sep='\t',dtype={'gene_id':'string'},low_memory=False)
 if len(sig)!=1326 or sig.gene_id.isna().any():raise ValueError('Frozen signature missing or unexpected size')
 sid=sig.gene_id.str.replace(r'\.\d+$','',regex=True)
 if sid.duplicated().any():raise ValueError('Signature IDs collapse to duplicates after version stripping')
 id_counts=ids.value_counts(dropna=True)
 id_to_row=pd.Series(np.arange(len(m)),index=ids).loc[~ids.duplicated(keep=False).to_numpy()]
 statuses=[]
 for i,base in enumerate(sid):
  n=int(id_counts.get(base,0))
  if n==0: status='absent'
  elif n>1: status='ambiguous_duplicate_id'
  else:
   idx=int(id_to_row.loc[base]);status='unique_complete' if eligible[idx] else ('unique_missing_counts' if anybad[idx] else 'unique_other_exclusion')
  statuses.append(status)
 mapping=sig[['gene_id']].copy();mapping['unversioned_id']=sid;mapping['mapping_status']=statuses
 if 'direction' in sig.columns:mapping['direction']=sig['direction']
 audit=m[['Ensembl.ID','Gene.symbol','Gene.ID']].copy()
 audit.insert(0,'row_index_zero_based',np.arange(len(m)))
 audit['unversioned_id']=ids
 audit['missing_or_invalid_count_samples']=invalid.sum(axis=1)
 audit['all_18_counts_invalid']=allbad
 audit['duplicate_unversioned_id']=duplicate
 audit['standard_ensembl_id']=standard
 audit['eligible_complete_unique_gene']=eligible
 qc=pd.DataFrame({'sample':SAMPLES,'missing_or_invalid_cells':invalid.sum(axis=0),'library_size_eligible':vals.loc[eligible,SAMPLES].sum(axis=0).astype('int64'),'detected_genes_eligible':(vals.loc[eligible,SAMPLES]>0).sum(axis=0)})
 OUT.mkdir(parents=True,exist_ok=False)
 audit.to_csv(OUT/'MATRIX_ROW_AUDIT_v4.tsv',sep='\t',index=False)
 mapping.to_csv(OUT/'SIGNATURE_MAPPING_v4.tsv',sep='\t',index=False)
 qc.to_csv(OUT/'COUNT_QC_v4.tsv',sep='\t',index=False)
 summary={'count_sha256':digest,'signature_sha256':sha(SIG),'matrix_rows':len(m),'invalid_any_rows':int(anybad.sum()),'invalid_all_18_rows':int(allbad.sum()),'partially_invalid_rows':int(partial.sum()),'duplicate_id_rows':int(duplicate.sum()),'eligible_unique_complete_rows':int(eligible.sum()),'signature_mapping':mapping.mapping_status.value_counts().to_dict(),'no_source_modification':True}
 (OUT/'AUDIT_MANIFEST_v4.json').write_text(json.dumps(summary,indent=2)+'\n')
 print('=== MATRIX AUDIT ===');print(json.dumps({k:v for k,v in summary.items() if k!='signature_mapping'},indent=2))
 print('=== SIGNATURE MAPPING ===');print(mapping.mapping_status.value_counts().to_string())
 if 'direction' in mapping:print(pd.crosstab(mapping.direction,mapping.mapping_status).to_string())
 print('=== COUNT QC ===');print(qc.to_string(index=False))
 print('Saved:',OUT)
if __name__=='__main__':main()
