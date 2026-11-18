#!/usr/bin/env python3
"""Make three provenance-aware evidence figures from TSV inputs.
No dataset/peak/network calls are invented; blanks are explicitly pending.
Requires pandas, matplotlib, numpy, networkx (only for PPI)."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CORE=['FDFT1','SQLE','CYP51A1','NSDHL','DHCR7']
OCC=['SQLE','FDFT1','CYP51A1','DHCR24','LSS']
C='#16375a'; BLUE='#4776a8'; ORANGE='#ba722b'

def save(fig,dest):
    dest.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(dest.with_suffix('.png'),dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(dest.with_suffix('.svg'),bbox_inches='tight',facecolor='white')
    plt.close(fig)

def fig1(data,out):
    d=pd.read_csv(data/'sterol_core.tsv',sep='\t')
    val=data/'external_validation.tsv'
    x=pd.read_csv(val,sep='\t') if val.exists() else pd.DataFrame(columns=['dataset','contrast','gene','log2FC','padj'])
    if not x.empty:
        x['log2FC']=pd.to_numeric(x.log2FC,errors='coerce')
        if x.duplicated(['dataset','contrast','gene']).any():raise ValueError('Duplicate external gene/contrast combinations')
    contexts=[('GSE267018','HeLa KO vs WT\n(primary)'),('GSE271000','HepG2 shSREBF2\n(external)'),('GSE147770','HAP1 KO vs WT\n(external)')]
    names=list(d.gene)
    arr=np.full((len(names),3),np.nan);arr[:,0]=d.primary_log2FC.values
    if not x.empty:
      for j,(ds,_) in enumerate(contexts[1:],1):
        choices=x[x.dataset==ds]
        if choices.contrast.nunique()>1:raise ValueError(f'{ds}: multiple contrasts. Filter to a prespecified contrast before plotting')
        m=dict(zip(choices.gene,choices.log2FC))
        arr[:,j]=[m.get(g,np.nan) for g in names]
    fig,ax=plt.subplots(figsize=(10.8,6.3))
    from matplotlib.colors import TwoSlopeNorm
    norm=TwoSlopeNorm(vcenter=0,vmin=-6,vmax=6)
    im=ax.imshow(np.ma.masked_invalid(arr),cmap='RdBu_r',norm=norm,aspect='auto')
    ax.set_facecolor('#eceff2')
    ax.set_yticks(range(len(names)),names)
    ax.set_xticks(range(3),[l for _,l in contexts],fontsize=10)
    ax.tick_params(length=0,pad=9)
    for i in range(len(names)):
      for j in range(3):
        v=arr[i,j]
        ax.text(j,i, f'{v:+.2f}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=9,
                color='white' if np.isfinite(v) and abs(v)>3.2 else '#263b52')
    ax.set_title('Figure A | 12-gene SREBF2 perturbation signature\nExternal effect estimates require matched, pre-specified contrasts',weight='bold',color=C,pad=15)
    cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.03);cb.set_label('log2 fold change (perturbation vs control)')
    fig.text(.5,.017,'NA = external gene-level DE results not provided; blank cells are NOT evidence of no effect.',ha='center',fontsize=9,color='#5e6b77')
    save(fig,out/'Figure_A_External_Validation')

def fig2(data,out):
    import networkx as nx
    edge_file=data/'ppi_edges.tsv'
    edges=pd.read_csv(edge_file,sep='\t')
    if not edges.empty:
      edges['score']=pd.to_numeric(edges.score,errors='raise')
      if not edges.score.between(0,1).all():raise ValueError('STRING score must be 0 to 1')
      if not edges.preferredName_A.isin(CORE).all() or not edges.preferredName_B.isin(CORE).all():
          raise ValueError('PPI edges must be limited to the specified 5 genes')
    G=nx.Graph();G.add_nodes_from(CORE)
    for _,e in edges.iterrows():G.add_edge(e.preferredName_A,e.preferredName_B,weight=float(e.score))
    pos={'FDFT1':(0,1),'SQLE':(1,1.4),'CYP51A1':(2,1),'NSDHL':(2,.0),'DHCR7':(0,.0)}
    fig,ax=plt.subplots(figsize=(9.8,5.6)); ax.axis('off')
    if G.number_of_edges():
        es=list(G.edges(data=True))
        nx.draw_networkx_edges(G,pos,ax=ax,edgelist=[(u,v) for u,v,_ in es],width=[1+5*d['weight'] for _,_,d in es],edge_color=BLUE,alpha=.75)
        nx.draw_networkx_edge_labels(G,pos,ax=ax,edge_labels={(u,v):f"{d['weight']:.2f}" for u,v,d in es},font_size=8)
    nx.draw_networkx_nodes(G,pos,ax=ax,node_color='#dbeaf6',edgecolors=C,node_size=3400,linewidths=1.4)
    nx.draw_networkx_labels(G,pos,ax=ax,font_size=12,font_weight='bold',font_color=C)
    ax.set_xlim(-.7,2.7);ax.set_ylim(-.55,1.9)
    ax.set_title('Figure B | Five-gene STRING association network',weight='bold',color=C,pad=15)
    msg=(f'{G.number_of_edges()} STRING associations displayed; edge labels are confidence scores (not physical binding).'
         if G.number_of_edges() else 'STRING edges NOT retrieved: node list only, no associations asserted.')
    fig.text(.5,.045,msg,ha='center',fontsize=10,color='#6d5363' if not G.number_of_edges() else C)
    save(fig,out/'Figure_B_STRING_Network')

def fig3(data,out):
    d=pd.read_csv(data/'occupancy_targets.tsv',sep='\t')
    cols=['GSE282800_HeLa_SREBP2','ENCSR000EZO_HepG2_SREBF2','GSE267019_HeLa_H3K27ac']
    if d.gene.tolist()!=OCC:raise ValueError('Unexpected target genes or ordering')
    colors={'not_tested':'#e1e8ef','overlap':'#2d6d63','no_overlap':'#d6a079'}
    fig,ax=plt.subplots(figsize=(11.2,4.8))
    for i,g in enumerate(OCC):
      for j,c in enumerate(cols):
        status=str(d.loc[i,c]);
        if status not in colors:raise ValueError('Invalid status: '+status)
        ax.add_patch(plt.Rectangle((j-.46,i-.41),.92,.82,fc=colors[status],ec='white',lw=1.2))
        ax.text(j,i,{'not_tested':'Not tested','overlap':'Peak overlap','no_overlap':'No peak overlap'}[status],ha='center',va='center',fontsize=9,color='white' if status=='overlap' else C)
    ax.set_xlim(-.5,2.5);ax.set_ylim(4.5,-.5)
    ax.set_yticks(range(5),OCC,fontsize=11)
    ax.set_xticks(range(3),['HeLa SREBP2\nGSE282800','HepG2 SREBF2\nENCODE ENCSR000EZO','HeLa H3K27ac\nGSE267019'],fontsize=10)
    ax.tick_params(length=0,pad=9)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_title('Figure C | Candidate regulatory-occupancy evidence audit\nPeak-level gene annotation pending',weight='bold',color=C,pad=17)
    fig.text(.5,.02,'Assay availability ≠ gene-specific occupancy; H3K27ac is not SREBP2 binding.',ha='center',fontsize=10,color='#5e6b77')
    save(fig,out/'Figure_C_Binding_Evidence_Audit')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=Path('data'));p.add_argument('--out',type=Path,default=Path('figures'))
    a=p.parse_args();fig1(a.data,a.out);fig2(a.data,a.out);fig3(a.data,a.out)
    print('Created 3 figures in PNG/SVG; pending evidence shown explicitly.')
