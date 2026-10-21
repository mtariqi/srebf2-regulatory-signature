#!/usr/bin/env python3
"""Phase 13D: reference chromatin-state enrichment tests."""

import csv
from collections import Counter
from pathlib import Path

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

INPUT = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)
OUTPUT = Path(
    "results/phase13/SREBF2_E117_STATE_FISHER_ENRICHMENT.tsv"
)

GROUPS = {
    "decreased_in_KO": 789,
    "increased_in_KO": 527,
    "background": 15781,
}

CONTRASTS = [
    ("decreased_in_KO", "background"),
    ("increased_in_KO", "background"),
    ("increased_in_KO", "decreased_in_KO"),
]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

counts = Counter()
totals = Counter()
seen = set()

with INPUT.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        gene = row["gene_id"]
        group = row["signature_group"]
        state = row["chromhmm_state"]

        if gene in seen:
            raise ValueError(f"Duplicate gene: {gene}")
        seen.add(gene)

        if group not in GROUPS:
            raise ValueError(f"Unknown group: {group}")

        totals[group] += 1
        counts[(group, state)] += 1

if dict(totals) != GROUPS:
    raise ValueError(f"Unexpected group counts: {totals}")

states = sorted(
    {
        state for group, state in counts
        if state != "UNANNOTATED"
    },
    key=lambda name: int(name.split("_", 1)[0]),
)

if len(states) != 15:
    raise ValueError("Expected exactly 15 chromatin states")

assigned = {
    group: totals[group] - counts[(group, "UNANNOTATED")]
    for group in GROUPS
}

results = []

for group_a, group_b in CONTRASTS:
    for state in states:
        a = counts[(group_a, state)]
        b = assigned[group_a] - a
        c = counts[(group_b, state)]
        d = assigned[group_b] - c

        odds_ratio, pvalue = fisher_exact(
            [[a, b], [c, d]],
            alternative="two-sided",
        )

        results.append({
            "contrast": f"{group_a}_vs_{group_b}",
            "state": state,
            "group_a_state": a,
            "group_a_other": b,
            "group_b_state": c,
            "group_b_other": d,
            "group_a_percent": 100 * a / assigned[group_a],
            "group_b_percent": 100 * c / assigned[group_b],
            "odds_ratio": odds_ratio,
            "pvalue": pvalue,
        })

adjusted = multipletests(
    [row["pvalue"] for row in results],
    method="fdr_bh",
)[1]

for row, qvalue in zip(results, adjusted):
    row["fdr_bh_45tests"] = qvalue

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("x", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(results[0]),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(results)

print("=== Phase 13D Fisher enrichment ===")
print("Status: PASS")
print("Tests:", len(results))
print("Assigned genes:", sum(assigned.values()))
print("Unannotated genes:", sum(
    counts[(g, "UNANNOTATED")] for g in GROUPS
))
print("FDR < 0.05:", sum(
    row["fdr_bh_45tests"] < 0.05 for row in results
))
print("Output:", OUTPUT)
