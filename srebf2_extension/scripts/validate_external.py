#!/usr/bin/env python3
"""QC and concordance summary for one contrast per external RNA-seq dataset.
Input table from re-run DESeq2: dataset,contrast,gene,log2FC,padj.
Does not substitute statin-treated cells for SREBF2 genetic perturbation."""
import argparse
from pathlib import Path
import numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--primary',default='data/sterol_core.tsv');p.add_argument('--external',required=True);p.add_argument('--out',default='data/external_validation_summary.tsv');a=p.parse_args()
a0=pd.read_csv(a.primary,sep='\t');ext=pd.read_csv(a.external,sep='\t');assert {'dataset','contrast','gene','log2FC','padj'}.issubset(ext)
if ext.empty:raise SystemExit('No external DE contrasts provided. No validation metrics can be reported.')
if ext.duplicated(['dataset','contrast','gene']).any():raise ValueError('Duplicate gene per dataset contrast')
for k in ['log2FC','padj']:ext[k]=pd.to_numeric(ext[k],errors='raise')
rows=[]
for (ds,contrast),sub in ext.groupby(['dataset','contrast']):
 m=a0.merge(sub,on='gene',how='left');valid=m.log2FC.notna();down=(m.loc[valid,'log2FC']<0)
 rows.append({'dataset':ds,'contrast':contrast,'n_12_tested':int(valid.sum()),'n_same_down_direction':int(down.sum()),'n_FDR05_same_direction':int(((m.loc[valid,'padj']<.05)&down).sum()),'median_external_log2FC':float(m.loc[valid,'log2FC'].median()) if valid.any() else np.nan,'spearman_primary_external':float(m.loc[valid,['primary_log2FC','log2FC']].corr(method='spearman').iloc[0,1]) if valid.sum()>=3 else np.nan})
out=pd.DataFrame(rows);Path(a.out).parent.mkdir(parents=True,exist_ok=True);out.to_csv(a.out,sep='\t',index=False);print(out.to_string(index=False))
