#!/usr/bin/env python3
"""Download human STRING association network using 12 SREBF2 sterol-core seeds + neighbors."""
import argparse,datetime,json
from pathlib import Path
import pandas as pd
import requests
SEEDS=['TM7SF2','CYP51A1','MSMO1','SQLE','FDFT1','HSD17B7','LSS','NSDHL','DHCR7','EBP','SC5D','DHCR24']
def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='srebf2_pca_gnn/data');p.add_argument('--add-nodes',type=int,default=70);p.add_argument('--score',type=int,default=700);args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    if args.add_nodes<30:raise SystemExit('Use at least 30 neighbor nodes for exploratory graph learning')
    url='https://string-db.org/api/tsv/network'
    params={'identifiers':'\r'.join(SEEDS),'species':9606,'required_score':args.score,'add_nodes':args.add_nodes,'network_type':'functional','caller_identity':'srebf2_regulatory_signature'}
    r=requests.post(url,data=params,timeout=90);r.raise_for_status()
    (out/'string_expanded_raw.tsv').write_text(r.text)
    df=pd.read_csv(out/'string_expanded_raw.tsv',sep='\t')
    expected=['preferredName_A','preferredName_B','score']
    if not all(x in df.columns for x in expected):raise ValueError('Unexpected STRING schema: '+str(df.columns.tolist()))
    df[expected].to_csv(out/'string_expanded_edges.tsv',sep='\t',index=False)
    meta={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':url,'parameters':params,'seed_genes':SEEDS,'returned_rows':len(df),'returned_nodes':len(set(df.preferredName_A)|set(df.preferredName_B))}
    (out/'string_expanded_provenance.json').write_text(json.dumps(meta,indent=2))
    print('PASS STRING expansion:',meta['returned_nodes'],'nodes,',meta['returned_rows'],'edges')
if __name__=='__main__':main()
