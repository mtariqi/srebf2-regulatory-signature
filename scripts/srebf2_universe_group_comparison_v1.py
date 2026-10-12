#!/usr/bin/env python3
"""Exploratory universe comparisons; use after promoter extraction."""
import csv
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import mannwhitneyu

BASE = Path('results/binding/GSE282800')
INPUT = BASE / 'SREBF2_DISCOVERY_UNIVERSE_PROMOTER_SIGNAL_v1.tsv'
OUT = BASE / 'DISCOVERY_UNIVERSE_GROUP_COMPARISONS_v1.tsv'


def median(x):
    return float(np.median(x)) if len(x) else float('nan')


def compare(a, b):
    x, y = np.array(a), np.array(b)
    u, p = mannwhitneyu(x, y, alternative='two-sided', method='asymptotic')
    return {'n_a': len(x), 'n_b': len(y), 'median_a': median(x),
            'median_b': median(y), 'median_difference': median(x) - median(y),
            'U': float(u), 'p_exploratory': float(p),
            'AUC': float(u / (len(x) * len(y)))}


with INPUT.open(newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
if len(rows) != 17158:
    raise ValueError(f'Expected 17158 rows, got {len(rows)}')

eligible = [r for r in rows if r['mapping_status'] == 'mapped_with_signal'
            and r['primary_nuclear_chromosome'] == 'True']
print('All evaluated genes:', len(rows))
print('Eligible nuclear genes:', len(eligible))
print('Excluded by status/chromosome:', len(rows) - len(eligible))
print('Eligible categories:', dict(Counter(
    r['direction'] if r['signature_selected'].lower() == 'true' else 'unselected'
    for r in eligible)))

for r in eligible:
    r['_group'] = (r['direction'] if r['signature_selected'].lower() == 'true'
                   else 'unselected')

metrics = ['FBS_mean', 'depletion_mean', 'depletion_minus_FBS']
comparisons = [('decreased_in_KO', 'unselected'),
               ('increased_in_KO', 'unselected'),
               ('decreased_in_KO', 'increased_in_KO')]
results = []
for metric in metrics:
    for a, b in comparisons:
        x = [float(r[metric]) for r in eligible if r['_group'] == a]
        y = [float(r[metric]) for r in eligible if r['_group'] == b]
        if not x or not y:
            raise ValueError(f'Empty group for {metric}: {a} vs {b}')
        result = {'metric': metric, 'group_a': a, 'group_b': b, **compare(x, y)}
        results.append(result)
        print(f'{metric}: {a} vs {b}; n={len(x)}/{len(y)}; '
              f'medians={result["median_a"]:.6f}/{result["median_b"]:.6f}; '
              f'AUC={result["AUC"]:.4f}; p={result["p_exploratory"]:.3g}')

with OUT.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(results)
print('Saved:', OUT)
print('CAUTION: exploratory unadjusted gene-level tests; not differential-binding tests.')
