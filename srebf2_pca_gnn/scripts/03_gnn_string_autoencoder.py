#!/usr/bin/env python3
"""Exploratory GCN graph autoencoder predicting *withheld STRING associations*.

Uses only training-graph topology + a frozen seed-gene indicator. No claim of
SREBF2 binding, expression prediction, or biological validation is made.
"""
import argparse,json,random
from pathlib import Path
import numpy as np,pandas as pd,networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch
from torch import nn
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.decomposition import PCA

SEEDS={'TM7SF2','CYP51A1','MSMO1','SQLE','FDFT1','HSD17B7','LSS','NSDHL','DHCR7','EBP','SC5D','DHCR24'}

def build_graph(path):
    df=pd.read_csv(path,sep='\t')
    req={'preferredName_A','preferredName_B','score'}
    if not req.issubset(df): raise ValueError('Expected columns: '+str(req))
    G=nx.Graph()
    for a,b,s in df[['preferredName_A','preferredName_B','score']].itertuples(index=False,name=None):
        if a==b:continue
        if not np.isfinite(float(s)) or float(s)<0:continue
        G.add_edge(str(a),str(b),score=float(s))
    if G.number_of_nodes()<30:raise ValueError('Graph too small (<30 nodes). Expand STRING network first; five-gene complete graph is invalid for held-out edge prediction.')
    n=G.number_of_nodes();possible=n*(n-1)//2
    if G.number_of_edges()>=possible or G.number_of_edges()<n:
        raise ValueError('Graph is complete or too sparse; cannot perform defensible evaluation.')
    return G

def split_edges(G,seed=13):
    rng=random.Random(seed)
    nodes=sorted(G.nodes());E={tuple(sorted(e)) for e in G.edges()}
    # Keep a spanning forest in the training network to preserve each connected component.
    required=set()
    for comp in nx.connected_components(G):
        for e in nx.minimum_spanning_tree(G.subgraph(comp)).edges(): required.add(tuple(sorted(e)))
    other=list(E-required);rng.shuffle(other)
    ntest=max(2,round(.15*len(E)));nval=max(2,round(.10*len(E)))
    if len(other)<ntest+nval:raise ValueError('Insufficient cycle edges for heldout split after preserving connectivity')
    test=other[:ntest];val=other[ntest:ntest+nval];train=list(E-set(test)-set(val))
    all_non=[(nodes[i],nodes[j]) for i in range(len(nodes)) for j in range(i+1,len(nodes)) if (nodes[i],nodes[j]) not in E]
    rng.shuffle(all_non)
    if len(all_non)<len(E):raise ValueError('Too few negative pairs: graph overly dense')
    ntest_neg=all_non[:len(test)];nval_neg=all_non[len(test):len(test)+len(val)]
    train_neg=all_non[len(test)+len(val):len(test)+len(val)+len(train)]
    if len(train_neg)<len(train):raise ValueError('Not enough train negatives')
    return nodes,train,test,val,train_neg,ntest_neg,nval_neg

class GCN(nn.Module):
    def __init__(self,d):
        super().__init__();self.l1=nn.Linear(d,32);self.l2=nn.Linear(32,16)
    def forward(self,x,adj):
        return self.l2(adj@torch.relu(self.l1(adj@x)))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--edges',default='srebf2_pca_gnn/data/string_expanded_edges.tsv')
    p.add_argument('--effects',default='srebf2_pca_gnn/data/hela_sterol_core_effects.tsv')
    p.add_argument('--out',default='srebf2_pca_gnn/figures')
    p.add_argument('--seed',type=int,default=13)
    p.add_argument('--epochs',type=int,default=160)
    args=p.parse_args()
    np.random.seed(args.seed);random.seed(args.seed);torch.manual_seed(args.seed)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    G=build_graph(args.edges)
    nodes,train,test,val,neg_train,neg_test,neg_val=split_edges(G,args.seed)
    idmap={s:i for i,s in enumerate(nodes)}
    T=nx.Graph();T.add_nodes_from(nodes);T.add_edges_from(train)
    n=len(nodes)
    A=nx.to_numpy_array(T,nodelist=nodes)
    A+=np.eye(n);deg=A.sum(axis=1)
    A=A/np.sqrt(deg[:,None]*deg[None,:])
    adj=torch.as_tensor(A,dtype=torch.float32)
    # Train-only structural features + prior status as one of the 12 curated cholesterol genes.
    degrees=np.array([T.degree(x) for x in nodes],float)
    clustering=np.array([nx.clustering(T,x) for x in nodes],float)
    seed=np.array([int(x in SEEDS) for x in nodes],float)
    feat=np.stack([np.log1p(degrees),clustering,seed,np.ones(n)],axis=1)
    mean=feat[:,:2].mean(0);std=feat[:,:2].std(0)+1e-8;feat[:,:2]=(feat[:,:2]-mean)/std
    x=torch.as_tensor(feat,dtype=torch.float32)
    train_idx=torch.tensor([[idmap[a],idmap[b]] for a,b in train+neg_train],dtype=torch.long)
    ytrain=torch.tensor([1.]*len(train)+[0.]*len(neg_train),dtype=torch.float32)
    model=GCN(feat.shape[1]);optimizer=torch.optim.Adam(model.parameters(),lr=.008,weight_decay=.001)
    for ep in range(args.epochs):
        model.train();optimizer.zero_grad();z=model(x,adj)
        preds=(z[train_idx[:,0]]*z[train_idx[:,1]]).sum(dim=1)
        loss=nn.functional.binary_cross_entropy_with_logits(preds,ytrain)
        loss.backward();optimizer.step()
    model.eval()
    with torch.no_grad(): z=model(x,adj).numpy()
    def probs(edges):
        logits=np.array([float(z[idmap[a]]@z[idmap[b]]) for a,b in edges]);return 1/(1+np.exp(-np.clip(logits,-50,50)))
    def metrics(pos,neg):
        y=np.array([1]*len(pos)+[0]*len(neg));scores=np.r_[probs(pos),probs(neg)]
        return {'roc_auc':float(roc_auc_score(y,scores)),'average_precision':float(average_precision_score(y,scores)),'n_positive':len(pos),'n_negative':len(neg)}
    # Common-neighbor topological baseline on EXACT SAME edge split.
    def cnmetrics(pos,neg):
        y=np.array([1]*len(pos)+[0]*len(neg))
        scores=np.array([len(list(nx.common_neighbors(T,a,b))) for a,b in pos+neg])
        return {'roc_auc':float(roc_auc_score(y,scores)),'average_precision':float(average_precision_score(y,scores))}
    report={'type':'exploratory_link_reconstruction','n_nodes':n,'n_edges_full':G.number_of_edges(),'n_training_edges':len(train),'split_seed':args.seed,'epochs':args.epochs,'model':{'validation':metrics(val,neg_val),'test':metrics(test,neg_test)},'common_neighbors_baseline':{'validation':cnmetrics(val,neg_val),'test':cnmetrics(test,neg_test)},'limitations':['STRING is an association database, not a direct physical-interaction gold standard','Negatives are non-edges, not verified biological non-interactions','Same network supplied both features and reconstruction targets; holdout evaluation only checks graph structure','Only one split tested; results should not be generalized','No RNA response labels were used for prediction','Prior pathway gene selection can inflate apparent biological coherence']}
    (out/'GNN_exploratory_metrics.json').write_text(json.dumps(report,indent=2))
    pd.DataFrame(z,index=nodes,columns=[f'embedding_{i+1}' for i in range(z.shape[1])]).rename_axis('gene').to_csv(out/'GNN_node_embeddings.tsv',sep='\t')
    effects=pd.read_csv(args.effects,sep='\t').set_index('gene')
    lfc=effects['HeLa_GSE267018_KO_vs_WT_log2FC'].to_dict()
    # PCA projection of learned embeddings; diagnostic view, NOT sample-level RNA-seq PCA.
    coords=PCA(n_components=2).fit_transform(z)
    fig,ax=plt.subplots(figsize=(9,7))
    keep=np.array([name in SEEDS for name in nodes]);
    ax.scatter(coords[~keep,0],coords[~keep,1],s=25,color='#bdc7d2',alpha=.6,label='STRING added neighbor')
    core_names=[nodes[i] for i in range(n) if keep[i]]
    core_idx=np.where(keep)[0]
    sc=ax.scatter(coords[core_idx,0],coords[core_idx,1],c=[lfc.get(nodes[i],np.nan) for i in core_idx],cmap='coolwarm',vmin=-5.5,vmax=0,s=135,edgecolor='#333',linewidth=.6,label='Sterol-core gene')
    for i in core_idx:ax.annotate(nodes[i],(coords[i,0],coords[i,1]),xytext=(4,4),textcoords='offset points',fontsize=9)
    fig.colorbar(sc,ax=ax,label='HeLa knockout vs WT log2FC (display overlay only)')
    ax.set_title('Exploratory STRING GCN embeddings',fontweight='bold')
    ax.set_xlabel('Embedding PC1');ax.set_ylabel('Embedding PC2')
    ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/'GNN_STRING_embedding.png',dpi=300);fig.savefig(out/'GNN_STRING_embedding.svg');plt.close(fig)
    pd.DataFrame({'gene':nodes,'embedding_PC1':coords[:,0],'embedding_PC2':coords[:,1],'is_sterol_core':keep,'hela_KO_log2FC_overlay':[lfc.get(q,np.nan) for q in nodes]}).to_csv(out/'GNN_embedding_projection.tsv',sep='\t',index=False)
    print(json.dumps(report,indent=2));print('Saved exploratory embeddings and diagnostic plot to',out)
if __name__=='__main__':main()
