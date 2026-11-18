#!/usr/bin/env python3
"""Strand-aware promoter interval overlaps for BED narrowPeak and hg38 gene TSS.
Input TSS TSV: gene,chrom,tss,strand (tss=1-based). Peak BED: 0-based half-open.
Does not assume peaks are matched to experimental KO and never infers causality."""
import argparse,gzip
from pathlib import Path
import pandas as pd
P=argparse.ArgumentParser();P.add_argument('--tss',required=True);P.add_argument('--peaks',required=True)
P.add_argument('--out',required=True);P.add_argument('--assay',required=True);P.add_argument('--upstream',type=int,default=2000);P.add_argument('--downstream',type=int,default=500)
a=P.parse_args();t=pd.read_csv(a.tss,sep='\t');assert {'gene','chrom','tss','strand'}.issubset(t)
assert a.upstream>=0 and a.downstream>=0
opens=gzip.open if a.peaks.endswith('.gz') else open
peaks={}
with opens(a.peaks,'rt') as h:
 for line in h:
  if line.startswith(('track','browser','#')):continue
  c=line.rstrip('\n').split('\t')
  if len(c)<3:continue
  try:s,e=int(c[1]),int(c[2])
  except ValueError:continue
  if s<0 or e<=s:raise ValueError('Invalid BED coordinates')
  peaks.setdefault(c[0],[]).append((s,e))
results=[]
for r in t.itertuples(index=False):
    if r.strand not in ['+','-']:raise ValueError('Strand must be + or -')
    p=int(r.tss)-1
    left,right=(p-a.upstream,p+a.downstream+1) if r.strand=='+' else (p-a.downstream,p+a.upstream+1)
    left=max(left,0)
    hit=[(s,e) for s,e in peaks.get(r.chrom,[]) if s<right and e>left]
    results.append({'gene':r.gene,'assay':a.assay,'chrom':r.chrom,'promoter_start0':left,'promoter_end0':right,'n_overlapping_peaks':len(hit),'peak_overlap':int(bool(hit))})
out=pd.DataFrame(results);Path(a.out).parent.mkdir(parents=True,exist_ok=True)
out.to_csv(a.out,sep='\t',index=False)
print('Wrote',len(out),'gene-level overlaps to',a.out)
