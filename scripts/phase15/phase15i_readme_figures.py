#!/usr/bin/env python3
"""Generate documented, publication-ready README figures from tracked TSVs."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

NAVY = '#183b65'; PALE = '#e9eef3'; INK = '#243449'; BLUE = '#6b8bac'

def plot_membership(root: Path):
    source=root/'results/phase15/figures/SREBF2_Figure5_Membership_Source.tsv'
    if not source.exists():
        raise FileNotFoundError(f'Missing tracked data: {source}')
    data=pd.read_csv(source,sep='\t')
    pathcols=[c for c in data.columns if c not in ('gene','number_of_pathways')]
    assert len(data)==12 and len(pathcols)==5 and data['gene'].is_unique
    vals=data[pathcols].to_numpy(dtype=int)
    assert np.isin(vals,[0,1]).all() and np.array_equal(vals.sum(axis=1),data['number_of_pathways'])
    order=data.sort_values(['number_of_pathways','gene'],ascending=[False,True]).reset_index(drop=True)
    vals=order[pathcols].to_numpy(dtype=int)
    labels=['GO\nsterol','GO\ncholesterol','Reactome\ncholesterol','Reactome\nSREBF','KEGG\nsteroid']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(12.2,8.1),facecolor='white')
    g=fig.add_gridspec(2,2,left=.18,right=.96,bottom=.18,top=.85,width_ratios=[5.4,1.1],height_ratios=[1.10,5],wspace=.12,hspace=.12)
    atop=fig.add_subplot(g[0,0]); main=fig.add_subplot(g[1,0]); side=fig.add_subplot(g[1,1],sharey=main)
    counts=vals.sum(axis=0)
    atop.bar(range(5),counts,width=.55,color=NAVY)
    atop.set_xlim(-.5,4.5); atop.set_ylim(0,14); atop.set_xticks([]); atop.set_ylabel('Core genes\nper pathway',color=INK)
    atop.spines[['top','right','bottom']].set_visible(False)
    for j,n in enumerate(counts): atop.text(j,n+.2,str(n),ha='center',weight='bold',color=INK)
    main.set_xlim(-.5,4.5); main.set_ylim(11.5,-.5)
    for i in range(12):
        for j in range(5):
            main.add_patch(Rectangle((j-.46,i-.41),.92,.82,color=PALE,linewidth=0))
            if vals[i,j]: main.add_patch(Rectangle((j-.32,i-.30),.64,.60,color=NAVY,linewidth=0))
    main.set_yticks(range(12),order['gene'],fontsize=11)
    main.set_xticks(range(5),labels,fontsize=10)
    main.tick_params(axis='both',length=0,pad=8); main.spines[:].set_visible(False)
    side.barh(range(12),vals.sum(axis=1),color=BLUE,height=.60)
    side.set_xlim(0,6.2); side.tick_params(axis='y',left=False,labelleft=False)
    side.set_xticks([]); side.spines[:].set_visible(False)
    for i,n in enumerate(vals.sum(axis=1)):side.text(n+.12,i,str(n),va='center',weight='bold',color=INK)
    side.set_title('Pathways\nper gene',color=INK,fontsize=10)
    fig.suptitle('Leading-edge membership of the 12-gene sterol core',y=.96,fontsize=16,weight='bold',color=INK)
    fig.legend(handles=[Patch(facecolor=NAVY,label='Member'),Patch(facecolor=PALE,label='Not a member')],loc='lower center',bbox_to_anchor=(.56,.085),ncol=2,frameon=False)
    fig.text(.56,.034,'Counts refer to these 12 genes only; GO, Reactome and KEGG annotations overlap.',ha='center',fontsize=9,color='#657487')
    target=root/'results/phase15/figures/SREBF2_Figure5_LeadingEdge_Membership'
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target.with_suffix('.png'),dpi=300,facecolor='white')
    fig.savefig(target.with_suffix('.svg'),facecolor='white')
    plt.close(fig)
    print(f'Figure 5: {target}.png / .svg')


def plot_forest(root: Path):
    src=root/'results/phase14/SREBF2_PHASE14A6_CHROMOSOME_EXACT_ENRICHMENT.tsv'
    if not src.exists():raise FileNotFoundError(f'Missing tracked data: {src}')
    d=pd.read_csv(src,sep='\t')
    expected={('increased_in_KO_vs_background','15_Quies'),('decreased_in_KO_vs_background','14_ReprPCWk')}
    assert expected.issubset(set(zip(d['comparison'],d['state'])))
    rows=d[d.apply(lambda r:(r['comparison'],r['state']) in expected,axis=1)].copy()
    rows=rows.sort_values('state',ascending=False).reset_index(drop=True)
    titles=[]
    for _,r in rows.iterrows():
        titles.append(('Increased in KO · quiescent E117 reference' if r['state']=='15_Quies' else 'Decreased in KO · weak Polycomb E117 reference'))
    fig,ax=plt.subplots(figsize=(10.8,3.7),facecolor='white')
    y=np.arange(len(rows))[::-1]
    for j,(_,r) in enumerate(rows.iterrows()):
        o=float(r['matched_or']);lo=float(r['ci95_lower']);hi=float(r['ci95_upper']); yy=y[j]
        ax.errorbar(o,yy,xerr=np.array([[o-lo],[hi-o]]),fmt='s',markersize=9,elinewidth=2.2,capsize=5,color=NAVY if r['state']=='15_Quies' else BLUE)
        ax.text(3.15,yy,f"OR {o:.2f} (95% CI {lo:.2f}–{hi:.2f})",va='center',fontsize=9,color=INK)
        ax.text(3.15,yy-.22,f"{int(r['n_pairs'])} pairs · BH q={float(r['fdr_bh_2tests']):.3g}",va='center',fontsize=8,color='#5e6d80')
    ax.axvline(1,ls='--',lw=1.3,color='#8997a8')
    ax.set_xlim(.4,9.0);ax.set_ylim(-.65,len(rows)-.30)
    ax.set_yticks(y,titles,fontsize=10);ax.set_xlabel('Matched odds ratio (log scale)')
    ax.set_xscale('log');ax.set_xticks([.5,1,2,4,8]);ax.set_xticklabels(['0.5','1','2','4','8'])
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=9)
    ax.set_title('Phase 14A6 · Chromosome-exact matched chromatin-state associations',fontsize=13,weight='bold',color=INK,pad=20)
    fig.subplots_adjust(left=.35,right=.98,top=.76,bottom=.19)
    target=root/'results/phase14/figures/SREBF2_Phase14A6_Matched_OR_Forest'
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target.with_suffix('.png'),dpi=300,facecolor='white')
    fig.savefig(target.with_suffix('.svg'),facecolor='white')
    plt.close(fig)
    print(f'Phase 14A6: {target}.png / .svg')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path.cwd())
    args=parser.parse_args()
    plot_membership(args.root)
    plot_forest(args.root)
