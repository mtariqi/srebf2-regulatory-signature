#!/usr/bin/env python3
"""Retrieve actual STRING v-current associations for five named Homo sapiens genes.
Requires internet on the user's machine; fails rather than fabricating edges.
Record retrieval date, API URL, and downloaded raw response alongside analysis."""
import argparse,datetime,io,json,urllib.parse,urllib.request
from pathlib import Path
import pandas as pd
GENES=['FDFT1','SQLE','CYP51A1','NSDHL','DHCR7']
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=Path('data/ppi_edges.tsv'))
parser.add_argument('--score',type=int,default=700);args=parser.parse_args()
query=urllib.parse.urlencode({'identifiers':'\r'.join(GENES),'species':9606,'required_score':args.score,'add_nodes':0,'network_type':'functional','caller_identity':'srebf2-regulatory-signature'})
url='https://string-db.org/api/tsv/network?'+query
try:
 with urllib.request.urlopen(url,timeout=35) as r:raw=r.read().decode('utf-8')
except Exception as exc:raise SystemExit('STRING retrieval unavailable; source table unchanged: '+str(exc))
df=pd.read_csv(io.StringIO(raw),sep='\t')
need=['preferredName_A','preferredName_B','score'];assert set(need).issubset(df)
df=df[df.preferredName_A.isin(GENES)&df.preferredName_B.isin(GENES)].copy()
if df.empty:raise SystemExit('No qualifying STRING edges returned; existing data unchanged.')
out=df[need].drop_duplicates().copy();out['source']='STRING API functional network; required_score='+str(args.score)
args.out.parent.mkdir(parents=True,exist_ok=True)
raw_path=args.out.with_suffix('.string_raw.tsv');raw_path.write_text(raw)
out.to_csv(args.out,sep='\t',index=False)
args.out.with_suffix('.provenance.json').write_text(json.dumps({'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':url,'taxon':9606,'min_score':args.score,'genes':GENES,'n_edges':len(out)},indent=2))
print('Saved verified STRING associations:',len(out),'edges to',args.out)
