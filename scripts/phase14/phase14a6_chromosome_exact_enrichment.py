#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase14/phase14a6_chromosome_exact_enrichment.py
# Purpose     : Phase 14A.6: paired enrichment after chromosome-exact matching
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 14A.6: paired enrichment after chromosome-exact matching."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

SOURCE = Path(
    "results/phase14/"
    "SREBF2_PHASE14A5_CHROMOSOME_EXACT_PAIRS.tsv"
)

ANNOTATION = Path(
    "results/phase13/"
    "SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)

OUTPUT = Path(
    "results/phase14/"
    "SREBF2_PHASE14A6_CHROMOSOME_EXACT_ENRICHMENT.tsv"
)

TESTS = [
    ("increased_in_KO_vs_background", "15_Quies", 425),
    ("decreased_in_KO_vs_background", "14_ReprPCWk", 677),
]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

pairs = pd.read_csv(SOURCE, sep="\t")
genes = pd.read_csv(
    ANNOTATION,
    sep="\t",
    usecols=["gene_id", "chrom", "chromhmm_state"],
)

if len(pairs) != 1102:
    raise ValueError("Unexpected pair count")

if not genes["gene_id"].is_unique:
    raise ValueError("Duplicate gene IDs")

lookup = genes.set_index("gene_id")

if pairs.duplicated(
    ["comparison", "signature_gene_id"]
).any():
    raise ValueError("Duplicate signature genes")

if pairs.duplicated(
    ["comparison", "background_gene_id"]
).any():
    raise ValueError("Background control reused")

if (
    pairs["standardized_distance"] > 0.25 + 1e-10
).any():
    raise ValueError("Matching caliper exceeded")

for prefix, id_column, state_column in [
    ("signature", "signature_gene_id", "signature_state"),
    ("background", "background_gene_id", "background_state"),
]:
    ids = pairs[id_column]
    observed_chrom = ids.map(lookup["chrom"])
    observed_state = ids.map(lookup["chromhmm_state"])

    if observed_chrom.isna().any() or observed_state.isna().any():
        raise ValueError(f"Missing {prefix} annotations")

    if not (
        observed_chrom.to_numpy()
        == pairs["chromosome"].to_numpy()
    ).all():
        raise ValueError(f"{prefix} chromosome mismatch")

    if not (
        observed_state.to_numpy()
        == pairs[state_column].to_numpy()
    ).all():
        raise ValueError(f"{prefix} state mismatch")

results = []

for comparison, state, expected_n in TESTS:
    sub = pairs.loc[
        pairs["comparison"] == comparison
    ].copy()

    if len(sub) != expected_n:
        raise ValueError(
            f"Unexpected pairs: {comparison}: {len(sub)}"
        )

    signature = sub["signature_state"].eq(state).to_numpy()
    background = sub["background_state"].eq(state).to_numpy()

    both = int(np.sum(signature & background))
    sig_only = int(np.sum(signature & ~background))
    bg_only = int(np.sum(~signature & background))
    neither = int(np.sum(~signature & ~background))

    if both + sig_only + bg_only + neither != len(sub):
        raise RuntimeError("Invalid contingency counts")

    discordant = sig_only + bg_only

    if discordant:
        test = binomtest(
            sig_only,
            discordant,
            p=0.5,
            alternative="two-sided",
        )
        pvalue = float(test.pvalue)
        interval = test.proportion_ci(
            confidence_level=0.95,
            method="exact",
        )

        def odds(probability):
            if probability == 0:
                return 0.0
            if probability == 1:
                return float("inf")
            return probability / (1 - probability)

        matched_or = (
            sig_only / bg_only
            if bg_only > 0
            else float("inf")
        )
        ci_lower = odds(interval.low)
        ci_upper = odds(interval.high)
    else:
        pvalue = 1.0
        matched_or = float("nan")
        ci_lower = float("nan")
        ci_upper = float("nan")

    results.append({
        "comparison": comparison,
        "state": state,
        "n_pairs": len(sub),
        "both_present": both,
        "signature_only": sig_only,
        "background_only": bg_only,
        "neither_present": neither,
        "discordant_pairs": discordant,
        "signature_state_percent": (
            100 * (both + sig_only) / len(sub)
        ),
        "background_state_percent": (
            100 * (both + bg_only) / len(sub)
        ),
        "matched_or": matched_or,
        "ci95_lower": ci_lower,
        "ci95_upper": ci_upper,
        "exact_pvalue": pvalue,
        "chromosomes_represented": sub["chromosome"].nunique(),
    })

out = pd.DataFrame(results)

out["fdr_bh_2tests"] = multipletests(
    out["exact_pvalue"],
    method="fdr_bh",
)[1]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

out.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    mode="x",
    lineterminator="\n",
)

print("=== Phase 14A.6 chromosome-exact enrichment ===")
print("Tests:", len(out))

for _, row in out.iterrows():
    print(f"\n=== {row['comparison']}: {row['state']} ===")
    print("Matched pairs:", row["n_pairs"])
    print("Chromosomes:", row["chromosomes_represented"])
    print("Both state-present:", row["both_present"])
    print("Signature-only:", row["signature_only"])
    print("Background-only:", row["background_only"])
    print("Neither state-present:", row["neither_present"])
    print("Discordant pairs:", row["discordant_pairs"])
    print("Matched OR:", f"{row['matched_or']:.4f}")
    print(
        "Exact 95% CI:",
        f"{row['ci95_lower']:.4f} - "
        f"{row['ci95_upper']:.4f}"
    )
    print("Exact p-value:", f"{row['exact_pvalue']:.6g}")
    print("BH FDR:", f"{row['fdr_bh_2tests']:.6g}")

print("\nPASS: Phase 14A.6 completed")
print("Output:", OUTPUT)
