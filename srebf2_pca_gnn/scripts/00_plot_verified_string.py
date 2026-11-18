#!/usr/bin/env python3
"""Plot empirical five-gene STRING network, without implying this is a GNN output."""
import argparse
from pathlib import Path
import pandas as pd,numpy as np,networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--edges',default='srebf2_pca_gnn/data/string_core_5_verified.tsv');p.add_argument('--effects',default='srebf2_pca_gnn/data/hela_sterol_core_effects.tsv');p.add_argument('--out',default='srebf2_pca_gnn/figures');a=p.parse_args()
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 d=pd.read_csv(a.edges,sep='\t');G=nx.Graph()
 for _,r in d.iterrows():G.add_edge(r.preferredName_A,r.preferredName_B,score=float(r.score))
 assert G.number_of_nodes()==5 and G.number_of_edges()==10, 'Unexpected verified core network shape'
 lfc=pd.read_csv(a.effects,sep='\t').set_index('gene')['HeLa_GSE267018_KO_vs_WT_log2FC'].to_dict()
 pos=nx.circular_layout(G,scale=1)
 fig,ax=plt.subplots(figsize=(8.8,7.3),facecolor='white')
 weights=[2+(G[u][v]['score']-.9)*30 for u,v in G.edges()]
 nx.draw_networkx_edges(G,pos,ax=ax,width=weights,alpha=.65,edge_color='#889daf')
 nodes=list(G.nodes())
 nx.draw_networkx_nodes(G,pos,ax=ax,nodelist=nodes,node_color=[lfc[n] for n in nodes],cmap='coolwarm',vmin=-5.5,vmax=0,node_size=2300,edgecolors='#23415a',linewidths=1.1)
 nx.draw_networkx_labels(G,pos,ax=ax,font_weight='bold',font_size=11,font_color='white')
 edge_labels={(u,v):f'{G[u][v]["score"]:.3f}' for u,v in G.edges()}
 nx.draw_networkx_edge_labels(G,pos,edge_labels,ax=ax,font_size=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':.3})
 ax.set_title('Measured STRING functional-association network',fontweight='bold',fontsize=15,pad=18)
 ax.text(.5,-.03,'5 curated sterol enzymes · 10 associations · score threshold ≥ 0.700\nNode shading = HeLa KO log2FC (descriptive overlay, not a graph-model result)',transform=ax.transAxes,ha='center',fontsize=9,color='#555')
 ax.axis('off');fig.tight_layout();fig.savefig(out/'STRING_5_gene_verified.png',dpi=300,bbox_inches='tight');fig.savefig(out/'STRING_5_gene_verified.svg',bbox_inches='tight');plt.close(fig)
 print('PASS: verified five-gene STRING diagram saved')
if __name__=='__main__':main()
