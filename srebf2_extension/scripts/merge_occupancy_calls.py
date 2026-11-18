#!/usr/bin/env python3
"""Update verified target-overlap statuses from one or more check_peak_overlap outputs.
Each --mapping is LABEL=FILE, where LABEL is a column in occupancy_targets.tsv.
Do not call 'no_overlap' for genes absent from the input gene list."""
import argparse
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--table',default='data/occupancy_targets.tsv');p.add_argument('--mapping',action='append',required=True);p.add_argument('--out',default='data/occupancy_targets_verified.tsv');a=p.parse_args()
d=pd.read_csv(a.table,sep='\t')
for spec in a.mapping:
 label,file=spec.split('=',1)
 if label not in d:raise ValueError('Unknown assay: '+label)
 t=pd.read_csv(file,sep='\t');assert {'gene','peak_overlap'}.issubset(t)
 if t.gene.duplicated().any():raise ValueError('Multiple rows per gene in '+file)
 if not t.peak_overlap.isin([0,1]).all():raise ValueError('peak_overlap must be binary')
 found=dict(zip(t.gene,t.peak_overlap))
 d[label]=[('overlap' if found[g]==1 else 'no_overlap') if g in found else old for g,old in zip(d.gene,d[label])]
out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);d.to_csv(out,sep='\t',index=False)
print('Verified occupancy calls saved:',out)
