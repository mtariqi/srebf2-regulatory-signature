#!/usr/bin/env python3
"""Sample-level SREBF2 HeLa RNA-seq PCA; logCPM is QC approximation, not DESeq2 VST."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def load_counts(path, sample_ids):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    if df.shape[1] <= 2:
        raise ValueError('Expected gene-by-sample count matrix; got fewer than 3 columns')
    # Detect exact runs or suffixes; require bijective match.
    mapping = {}
    for run in sample_ids:
        matches = [c for c in df.columns if c == run or str(c).split('/')[-1].split('.')[0] == run]
        if len(matches) != 1:
            raise ValueError(f'Missing/ambiguous count column {run}: {matches}')
        mapping[run] = matches[0]
    v = df[[mapping[x] for x in sample_ids]].apply(pd.to_numeric, errors='coerce')
    if v.isna().any().any() or (v < 0).any().any():
        raise ValueError('Counts have missing, nonnumeric or negative values')
    if (v.sum(axis=0) <= 0).any():
        raise ValueError('At least one sample has zero library size')
    v.columns = sample_ids
    return v


def compute_plot(counts, metadata, prefix, n_top, drop_run=None):
    md = metadata.copy()
    if drop_run:
        md = md.loc[md.run != drop_run].reset_index(drop=True)
    x = counts.loc[:, md.run].copy()
    lib = x.sum(axis=0)
    detected = (x >= 10).sum(axis=1) >= max(2, int(np.ceil(.25 * x.shape[1])))
    x = x.loc[detected]
    if x.shape[0] < 100:
        raise ValueError(f'Only {len(x)} expressed genes after filtering: check count matrix orientation')
    # Library size-normalized logCPM, intended only as exploratory QC PCA.
    logcpm = np.log2(x.div(lib, axis=1) * 1e6 + 1)
    var = logcpm.var(axis=1)
    top = logcpm.loc[var.nlargest(min(n_top, len(logcpm))).index]
    sample_matrix = top.T.to_numpy(dtype=float, copy=True)
    sample_matrix -= sample_matrix.mean(axis=0, keepdims=True)
    pca = PCA(n_components=min(3, sample_matrix.shape[0] - 1), svd_solver='full')
    pcs = pca.fit_transform(sample_matrix)
    scores = md.copy()
    for i in range(pcs.shape[1]): scores[f'PC{i+1}'] = pcs[:,i]
    scores.to_csv(str(prefix)+'_scores.tsv', sep='\t', index=False)
    pd.DataFrame({'PC':range(1,len(pca.explained_variance_ratio_)+1),'variance_fraction':pca.explained_variance_ratio_}).to_csv(str(prefix)+'_variance.tsv',sep='\t',index=False)
    fig, ax = plt.subplots(figsize=(9.1,6.6))
    colors = {'KO':'#ba4a40','WT':'#1d6484'}
    markers = {'FBS':'o','depletion':'s'}
    for _, row in scores.iterrows():
        ax.scatter(row.PC1, row.PC2, c=colors.get(row.genotype,'#777777'), marker=markers.get(row.treatment,'o'), s=90, alpha=.9, edgecolor='white',linewidth=.7)
        ax.annotate(str(row.run)[-4:], (row.PC1,row.PC2), xytext=(5,4),textcoords='offset points',fontsize=8)
    from matplotlib.lines import Line2D
    handles = [Line2D([],[], marker='o',linestyle='',label='KO',markerfacecolor=colors['KO'],markeredgecolor='none',markersize=9),Line2D([],[],marker='o',linestyle='',label='WT',markerfacecolor=colors['WT'],markeredgecolor='none',markersize=9),Line2D([],[],marker='o',linestyle='',color='#333',label='FBS',markersize=8),Line2D([],[],marker='s',linestyle='',color='#333',label='Cholesterol depletion',markersize=8)]
    ax.legend(handles=handles,frameon=False,loc='best')
    ax.axhline(0,color='#ddd',lw=.8);ax.axvline(0,color='#ddd',lw=.8)
    ax.set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% of variance)')
    ax.set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% of variance)')
    fig.suptitle('HeLa SREBF2 RNA-seq: sample-level PCA',fontsize=15,fontweight='bold')
    ax.set_title(f'GSE267018 | {len(scores)} samples | {min(n_top,len(logcpm)):,} variable expressed genes | logCPM',fontsize=10)
    fig.tight_layout(); fig.savefig(str(prefix)+'.png',dpi=300);fig.savefig(str(prefix)+'.svg');plt.close(fig)
    return len(scores), len(logcpm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts',default='results/counts/GSE267018_raw_count_matrix.tsv')
    ap.add_argument('--metadata',default='metadata/GSE267018_analysis_samples.tsv')
    ap.add_argument('--out',default='srebf2_pca_gnn/figures')
    ap.add_argument('--top-genes',type=int,default=2000)
    args=ap.parse_args()
    output=Path(args.out);output.mkdir(parents=True,exist_ok=True)
    metadata=pd.read_csv(args.metadata,sep='\t')
    required={'run','genotype','treatment'}
    if not required.issubset(metadata):raise ValueError('Metadata must have '+str(required))
    if metadata.run.duplicated().any():raise ValueError('Duplicate run ids in metadata')
    counts=load_counts(args.counts,metadata.run.tolist())
    for label,exclude in [('all12',None),('exclude_SRR28966285','SRR28966285')]:
        n,g=compute_plot(counts,metadata,output/f'PCA_{label}',args.top_genes,exclude)
        print(f'PASS: PCA_{label} from {n} actual samples and {g} expressed genes')
    print('Interpretation: exploratory sample QC only; logCPM is not variance stabilization and PCs are not tests of genotype effects.')
if __name__=='__main__':main()
