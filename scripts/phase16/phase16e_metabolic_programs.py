#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SREBF2 Regulatory Signature — Phase 16E
=======================================
Script       : phase16e_metabolic_programs.py
Project      : Historical chromatin-state context of transcriptional
               responses to SREBF2 knockout: an integrative reanalysis
               of HeLa RNA-seq
Repository   : https://github.com/mtariqi/srebf2-regulatory-signature
Author       : Md Tariqul Islam (project attribution; review before release)
Created      : 2026-11-19
Language     : Python 3
Dependencies : numpy, matplotlib (standard library: argparse, csv, pathlib)

Purpose
-------
Generate the manuscript Section 4.2 two-panel metabolic-programs figure:
(A) original negative GSEA NES and BH FDR comparisons; and
(B) a conceptual respiratory-chain gene-family schematic.

Data provenance and interpretation
----------------------------------
Default pathway statistics are transcribed from the V7–V9 manuscript and
are NOT recalculated from GSEA RDS objects. Pass --table with an audited
TSV to replace the frozen values. Panel B is descriptive, not an
expression heatmap or a verified count by respiratory-chain complex.
The 57-gene leading-edge total is manuscript-reported. The figure does not
establish oxygen consumption, ATP generation, sterol flux, causality, or
independent biological replication. Reactome translation terms are omitted
until their numerical results are verified.

Usage
-----
python3 scripts/phase16/phase16e_metabolic_programs.py \
  --output-prefix results/phase16/figures/Figure_Metabolic_Programs

Expected input TSV columns: pathway, NES, BH_FDR, source
Outputs: <prefix>.png, <prefix>.pdf, <prefix>_source.tsv
"""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

FROZEN = [
    ('Reactome · Cholesterol biosynthesis', -3.028, 6.06e-9, 'Reactome'),
    ('GO BP · Sterol biosynthetic process', -2.948, 4.05e-8, 'GO BP'),
    ('GO BP · Cholesterol biosynthetic process', -2.910, 4.05e-8, 'GO BP'),
    ('KEGG · Oxidative phosphorylation', -2.188, 1.04e-6, 'KEGG'),
]
FAMILIES = [('Complex I', 'NDUF'), ('Complex II', 'SDH'),
            ('Complex III', 'UQCR'), ('Complex IV', 'COX'),
            ('ATP synthase', 'ATP5')]

def read_rows(file):
    if file is None:
        return FROZEN
    with open(file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        expected = {'pathway', 'NES', 'BH_FDR', 'source'}
        if not expected.issubset(set(reader.fieldnames or [])):
            raise ValueError(f'TSV needs columns: {sorted(expected)}')
        rows = []
        for r in reader:
            nes, fdr = float(r['NES']), float(r['BH_FDR'])
            if not np.isfinite(nes) or not np.isfinite(fdr) or nes >= 0 or not (0 <= fdr <= 1):
                raise ValueError('Expected negative finite NES and finite BH_FDR in [0,1]')
            rows.append((r['pathway'], nes, fdr, r['source']))
    if not rows:
        raise ValueError('Input table is empty')
    return rows

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--table', type=Path, default=None, help='Audited pathway source TSV')
    p.add_argument('--output-prefix', type=Path,
                   default=Path('results/phase16/figures/Figure_Metabolic_Programs'))
    p.add_argument('--dpi', type=int, default=300)
    args = p.parse_args()
    data = read_rows(args.table)
    fig = plt.figure(figsize=(12.8, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.6, 1], hspace=.48)
    ax = fig.add_subplot(gs[0])
    names = [r[0] for r in data]
    vals = [-r[1] for r in data]
    ys = np.arange(len(data))
    ax.barh(ys, vals, height=.58, color='#277bb5')
    ax.set_yticks(ys, names, fontsize=10.5)
    ax.invert_yaxis()
    ax.set_xlabel('Magnitude of negative normalized enrichment score (−NES)')
    ax.set_title('A  Negative enrichment of sterol and respiratory pathways',
                 loc='left', fontsize=13, fontweight='bold', pad=18)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.grid(axis='x', alpha=.18)
    ax.set_axisbelow(True)
    for idx, (_, nes, fdr, _) in enumerate(data):
        ax.text(-nes + .07, idx, f'NES {nes:.3f}  BH FDR {fdr:.2e}', va='center', fontsize=9)
    ax.set_xlim(0, max(vals) + 2.65)

    ax2 = fig.add_subplot(gs[1])
    ax2.axis('off')
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 3.5)
    ax2.set_title('B  Respiratory-chain family schematic (descriptive)',
                  loc='left', fontsize=13, fontweight='bold', pad=12)
    for x, (complex_name, family) in zip(np.linspace(1.05, 8.95, 5), FAMILIES):
        ax2.add_patch(FancyBboxPatch((x-.86, 1.6), 1.72, 1.,
                         boxstyle='round,pad=0.12,rounding_size=0.12',
                         linewidth=1.2, edgecolor='#707070', facecolor='#f1f1f1'))
        ax2.text(x, 2.10, f'{complex_name}\n{family}', ha='center', va='center', fontsize=10)
    ax2.text(5, 1.11, '57 oxidative-phosphorylation leading-edge genes (manuscript-reported)',
             ha='center', va='center', fontsize=10.5, fontweight='bold')
    ax2.text(5, .43,
             'No measured respiratory activity, oxygen consumption, ATP generation or sterol flux.\n'
             'Reactome translation terms excluded: exact pathway statistics not supplied.',
             ha='center', va='center', fontsize=9.5, color='#444444')
    fig.suptitle('SREBF2 knockout in HeLa: broader metabolic transcriptional response',
                 fontsize=16, fontweight='bold', y=.99)
    fig.text(.5, .025,
             'Frozen manuscript GSEA estimates; overlapping annotations do not represent independent replication.',
             ha='center', fontsize=9, color='#555555')
    fig.subplots_adjust(left=.31, right=.97, top=.88, bottom=.09)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'pdf']:
        filename = args.output_prefix.with_suffix('.' + ext)
        fig.savefig(filename, dpi=args.dpi, bbox_inches='tight')
        print('CREATED', filename)
    plt.close(fig)
    source = args.output_prefix.with_name(args.output_prefix.name + '_source.tsv')
    with source.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pathway', 'NES', 'BH_FDR', 'source'])
        writer.writerows(data)
    print('CREATED', source)

if __name__ == '__main__':
    main()
