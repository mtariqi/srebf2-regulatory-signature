#!/usr/bin/env python3
"""Phase 14A.2: exact paired enrichment after background matching."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

SOURCE = Path(
    "results/phase14/SREBF2_PHASE14A2_MATCHED_PAIRS.tsv"
)

OUTPUT = Path(
    "results/phase14/SREBF2_PHASE14A2_MATCHED_ENRICHMENT.tsv"
)

TESTS = [
    ("increased_in_KO_vs_background", "15_Quies"),
    ("decreased_in_KO_vs_background", "14_ReprPCWk"),
]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

df = pd.read_csv(SOURCE, sep="\t")

if len(df) != 1293:
    raise ValueError("Unexpected matched-pair count")

if df[
    ["comparison", "signature_gene_id"]
].duplicated().any():
    raise ValueError("Duplicate signature gene within comparison")

if df[
    ["comparison", "background_gene_id"]
].duplicated().any():
    raise ValueError("Control reused within comparison")

if df["standardized_distance"].gt(0.25 + 1e-10).any():
    raise ValueError("Matching caliper violation")

results = []

for comparison, state in TESTS:

    sub = df[df["comparison"] == comparison].copy()

    expected_n = (
        514 if comparison.startswith("increased") else 779
    )

    if len(sub) != expected_n:
        raise ValueError(
            f"Unexpected pair count for {comparison}"
        )

    signature = sub["signature_state"].eq(state).to_numpy()
    control = sub["background_state"].eq(state).to_numpy()

    both_present = int((signature & control).sum())
    signature_only = int((signature & ~control).sum())
    control_only = int((~signature & control).sum())
    neither = int((~signature & ~control).sum())

    discordant = signature_only + control_only

    if (
        both_present + signature_only +
        control_only + neither != len(sub)
    ):
        raise RuntimeError("Pair counts do not sum correctly")

    if discordant == 0:
        raise RuntimeError(
            f"No discordant pairs for {comparison}, {state}"
        )

    # Exact McNemar p-value via the conditional binomial test.
    test = binomtest(
        signature_only,
        discordant,
        p=0.5,
        alternative="two-sided",
    )

    # Exact binomial CI, transformed to the matched odds-ratio scale.
    ci = test.proportion_ci(
        confidence_level=0.95,
        method="exact",
    )

    def odds(p):
        if p <= 0:
            return 0.0
        if p >= 1:
            return float("inf")
        return p / (1 - p)

    matched_or = (
        signature_only / control_only
        if control_only else float("inf")
    )

    results.append({
        "comparison": comparison,
        "state": state,
        "n_pairs": len(sub),
        "both_state_present": both_present,
        "signature_only": signature_only,
        "background_only": control_only,
        "neither_state_present": neither,
        "discordant_pairs": discordant,
        "signature_state_percent": (
            100 * (both_present + signature_only) / len(sub)
        ),
        "background_state_percent": (
            100 * (both_present + control_only) / len(sub)
        ),
        "matched_odds_ratio": matched_or,
        "exact_ci95_lower": odds(ci.low),
        "exact_ci95_upper": odds(ci.high),
        "pvalue_exact": test.pvalue,
    })

out = pd.DataFrame(results)

out["fdr_bh_2tests"] = multipletests(
    out["pvalue_exact"],
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

print("=== Phase 14A.2 matched enrichment ===")
print("Tests:", len(out))

for _, row in out.iterrows():
    print(f"\nComparison: {row['comparison']}")
    print("State:", row["state"])
    print("Matched pairs:", row["n_pairs"])
    print("Signature-only:", row["signature_only"])
    print("Background-only:", row["background_only"])
    print("Discordant pairs:", row["discordant_pairs"])

    print(
        "Matched OR:",
        f"{row['matched_odds_ratio']:.4f}"
    )

    print(
        "Exact 95% CI:",
        f"{row['exact_ci95_lower']:.4f}",
        "-",
        f"{row['exact_ci95_upper']:.4f}"
    )

    print(
        "Exact p-value:",
        f"{row['pvalue_exact']:.6g}"
    )

    print(
        "BH FDR (2 tests):",
        f"{row['fdr_bh_2tests']:.6g}"
    )

print("\nPASS: Paired enrichment calculations completed")
print("Output:", OUTPUT)
