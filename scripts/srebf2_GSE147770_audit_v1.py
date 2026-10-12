#!/usr/bin/env python3
"""Read-only Phase 6A GSE147770 metadata, count QC, and frozen signature mapping audit."""
import csv, gzip, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path.cwd()
COUNTS=ROOT/'data/raw/validation/GSE147770/GSE147770_STAR_readCounts.txt.gz'
MATRIX=ROOT/'metadata/validation/GSE147770/GSE147770_series_matrix.txt.gz'
SIGN=ROOT/'results/signature/discovery_v1/SREBF2_SIGNATURE_ALL.tsv'
OUT=ROOT/'results/validation/GSE147770/audit_v1'
SAMPLES=[f'{g}{i}{condition}' for g in 'CSW' for i in range(1,4) for condition in ('min','plus')]
EXPECTED_HASH='3f4a2307ede59793b77ac50c38fa519ece012c0a2128fe56ea4c0f44f0b7d598'

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    for p in (COUNTS,MATRIX,SIGN):
        if not p.is_file():raise FileNotFoundError(p)
    if OUT.exists() and any(OUT.iterdir()):raise FileExistsError(f'Refusing to overwrite: {OUT}')
    if digest(COUNTS)!=EXPECTED_HASH:raise ValueError('Count matrix checksum differs from downloaded matrix; audit source before proceeding')
    with gzip.open(MATRIX,'rt') as f:
        records={}
        for line in f:
            if line.startswith(('!Sample_title\t','!Sample_geo_accession\t','!Sample_characteristics_ch1\t')):
                fields=next(csv.reader([line.rstrip('\n')],delimiter='\t',quotechar='"'))
                records.setdefault(fields[0],[]).append(fields[1:])
    titles=records.get('!Sample_title',[[]])[0]
    gsm=records.get('!Sample_geo_accession',[[]])[0]
    if len(titles)!=18 or len(gsm)!=18:raise ValueError('Expected 18 GEO sample titles and accessions')
    sample_rows=[]
    for title,accession in zip(titles,gsm):
        m=re.match(r'^(?P<sample>[CSW][123](?:plus|min)):',title)
        if not m:raise ValueError(f'Unrecognized sample title: {title}')
        sample=m.group('sample'); group=sample[0]
        genotype={'W':'WT','S':'SREBF2_KO','C':'C12orf49_KO'}[group]
        condition='plus_FBS' if sample.endswith('plus') else 'minus_FBS'
        expected_label={'W':'WT','S':'SREBF2','C':'C12orf49'}[group]
        if expected_label.lower() not in title.lower():raise ValueError(f'Title/genotype mismatch: {title}')
        sample_rows.append(dict(sample=sample,gsm=accession,genotype=genotype,condition=condition,replicate=int(sample[1]),geo_title=title,primary_comparison=group in 'SW'))
    metadata=pd.DataFrame(sample_rows)
    if len(set(metadata['sample']))!=18:raise ValueError('Duplicate sample title IDs')
    expected_counts=metadata.groupby(['genotype','condition']).size()
    if not (expected_counts==3).all() or len(expected_counts)!=6:raise ValueError('Incorrect sample balance')
    with gzip.open(COUNTS,'rt') as f: header=f.readline().rstrip('\n').split('\t')
    if header[:3]!=['Ensembl.ID','Gene.symbol','Gene.ID']:raise ValueError(f'Unexpected matrix annotations: {header[:3]}')
    if len(header)!=21 or set(header[3:])!=set(SAMPLES):raise ValueError('Unexpected matrix sample columns')
    matrix=pd.read_csv(COUNTS,sep='\t',dtype={'Ensembl.ID':str,'Gene.symbol':str,'Gene.ID':str},low_memory=False)
    if len(matrix)==0:raise ValueError('Empty count matrix')
    ids=matrix['Ensembl.ID'].astype(str)
    if ids.isna().any() or ids.duplicated().any():raise ValueError('Missing or duplicate validation Ensembl IDs')
    if not ids.str.fullmatch(r'ENSG\d+(?:\.\d+)?').all():raise ValueError('Non-Ensembl IDs present')
    vals=matrix[header[3:]].apply(pd.to_numeric,errors='raise')
    arr=vals.to_numpy(dtype=float)
    if not np.isfinite(arr).all() or (arr<0).any() or not np.equal(arr,np.floor(arr)).all():raise ValueError('Nonfinite, negative or noninteger counts')
    qc=[]
    for sample in header[3:]:
        x=vals[sample]
        qc.append(dict(sample=sample,library_size=int(x.sum()),detected_genes=int((x>0).sum()),zero_fraction=float((x==0).mean())))
    signature=pd.read_csv(SIGN,sep='\t',dtype=str)
    if len(signature)!=1326 or 'gene_id' not in signature or 'direction' not in signature:raise ValueError('Unexpected frozen signature format')
    signature['ensembl_unversioned']=signature['gene_id'].str.replace(r'\.\d+$','',regex=True)
    if signature['ensembl_unversioned'].duplicated().any():raise ValueError('Ambiguous version-stripped signature IDs')
    map_count=ids.value_counts()
    signature['validation_row_count']=signature['ensembl_unversioned'].map(map_count).fillna(0).astype(int)
    signature['mapping_status']=np.select([signature.validation_row_count.eq(1),signature.validation_row_count.gt(1)],['unique','ambiguous'],default='absent')
    if signature.mapping_status.eq('ambiguous').any():raise ValueError('Ambiguous validation ID mapping')
    directions=signature.groupby(['direction','mapping_status']).size().reset_index(name='genes')
    manifest=dict(dataset='GSE147770',stage='phase_6A_audit_v1',counts_sha256=digest(COUNTS),geo_matrix_sha256=digest(MATRIX),signature_sha256=digest(SIGN),matrix_rows=len(matrix),samples=18,signature_total=1326,signature_unique_mapped=int(signature.mapping_status.eq('unique').sum()),signature_absent=int(signature.mapping_status.eq('absent').sum()),validation_identifier='unversioned Ensembl',mapping='remove terminal Ensembl version suffix only',no_analysis_performed=True)
    OUT.mkdir(parents=True,exist_ok=True)
    metadata.to_csv(OUT/'GSE147770_SAMPLE_METADATA_v1.tsv',sep='\t',index=False)
    pd.DataFrame(qc).merge(metadata[['sample','genotype','condition']],on='sample').to_csv(OUT/'GSE147770_COUNT_QC_v1.tsv',sep='\t',index=False)
    signature.to_csv(OUT/'GSE147770_GENE_MAPPING_AUDIT_v1.tsv',sep='\t',index=False)
    (OUT/'GSE147770_AUDIT_MANIFEST_v1.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('=== SAMPLE GROUPS ===')
    print(expected_counts.to_string())
    print('\n=== SIGNATURE MAPPING ===')
    print(directions.to_string(index=False))
    print('\n=== COUNT QC ===')
    print(pd.DataFrame(qc).to_string(index=False))
    print('\nSaved:',OUT)

if __name__=='__main__':main()
