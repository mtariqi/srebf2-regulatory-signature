#!/usr/bin/env python3
"""Read-only count-value diagnostic for GSE147770; never changes source or previous audits."""
from pathlib import Path
import gzip, hashlib, json
import numpy as np
import pandas as pd

ROOT=Path.cwd()
SRC=ROOT/'data/raw/validation/GSE147770/GSE147770_STAR_readCounts.txt.gz'
OUT=ROOT/'results/validation/GSE147770/count_diagnostic_v3'
EXPECTED='3f4a2307ede59793b77ac50c38fa519ece012c0a2128fe56ea4c0f44f0b7d598'

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def main():
    if sha256(SRC)!=EXPECTED:raise ValueError('Count file SHA256 mismatch')
    if OUT.exists() and any(OUT.iterdir()):raise FileExistsError(f'Refusing overwrite: {OUT}')
    with gzip.open(SRC,'rt') as f: header=f.readline().rstrip('\r\n').split('\t')
    if header[:3]!=['Ensembl.ID','Gene.symbol','Gene.ID'] or len(header)!=21:
        raise ValueError('Unexpected column layout')
    samples=header[3:]
    df=pd.read_csv(SRC,sep='\t',dtype=str,keep_default_na=False,na_filter=False)
    if list(df.columns)!=header:raise ValueError('Parsed columns differ from header')
    bad=[]; summary=[]
    for sample in samples:
        raw=df[sample].astype(str).str.strip()
        num=pd.to_numeric(raw.replace('',np.nan),errors='coerce')
        finite=np.isfinite(num.to_numpy(dtype=float))
        negative=finite & (num.to_numpy(dtype=float)<0)
        fractional=finite & (num.to_numpy(dtype=float)!=np.floor(num.to_numpy(dtype=float)))
        invalid=~finite | negative | fractional
        for ix in np.flatnonzero(invalid):
            reason=('missing_or_nonnumeric' if not finite[ix] else 'negative' if negative[ix] else 'fractional')
            bad.append({'row_index_zero_based':int(ix),'matrix_line_one_based':int(ix+2),
                        'Ensembl.ID':df.at[ix,'Ensembl.ID'],'Gene.symbol':df.at[ix,'Gene.symbol'],
                        'Gene.ID':df.at[ix,'Gene.ID'],'sample':sample,'raw_value':df.at[ix,sample],
                        'parsed_value':str(num.iloc[ix]),'reason':reason})
        summary.append({'sample':sample,'total_rows':len(df),'invalid_values':int(invalid.sum()),
                        'missing_or_nonnumeric':int((~finite).sum()),'negative':int(negative.sum()),
                        'fractional':int(fractional.sum())})
    OUT.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(summary).to_csv(OUT/'COUNT_VALUE_SUMMARY_v3.tsv',sep='\t',index=False)
    pd.DataFrame(bad,columns=['row_index_zero_based','matrix_line_one_based','Ensembl.ID','Gene.symbol','Gene.ID','sample','raw_value','parsed_value','reason']).to_csv(OUT/'COUNT_VALUE_EXCEPTIONS_v3.tsv',sep='\t',index=False)
    manifest={'dataset':'GSE147770','stage':'count_diagnostic_v3','counts_sha256':sha256(SRC),
              'matrix_rows':len(df),'samples':len(samples),'invalid_cells':len(bad),
              'source_modified':False,'repair_applied':False,
              'interpretation':'Raw count matrix must be assessed before DESeq2; no automatic rounding or imputation.'}
    (OUT/'COUNT_DIAGNOSTIC_MANIFEST_v3.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('=== COUNT VALUE SUMMARY ===')
    print(pd.DataFrame(summary).to_string(index=False))
    print('\n=== EXCEPTIONS (first 30) ===')
    print(pd.DataFrame(bad).head(30).to_string(index=False) if bad else 'None')
    print('\nTotal invalid count cells:',len(bad))
    print('Saved:',OUT)

if __name__=='__main__':main()
