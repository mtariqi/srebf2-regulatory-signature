#!/usr/bin/env python3
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase14/phase14a5_chromosome_exact_matching.py
# Purpose     : Phase 14A.5: chromosome- and biotype-exact matching
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================

"""Phase 14A.5: chromosome- and biotype-exact matching.

Deterministic greedy nearest-neighbor matching without replacement,
using a 0.25 standardized Euclidean distance caliper.

This is a sensitivity analysis, not globally optimal matching.
ChromHMM state is excluded from match selection.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SOURCE = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)
OUTPUT = Path(
    "results/phase14/"
    "SREBF2_PHASE14A5_CHROMOSOME_EXACT_PAIRS.tsv"
)

CALIPER = 0.25
FEATURES = ["log_expression", "log_length"]
GROUPS = ["increased_in_KO", "decreased_in_KO"]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

df = pd.read_csv(SOURCE, sep="\t")

if len(df) != 17097 or not df["gene_id"].is_unique:
    raise ValueError("Unexpected source gene universe")

df = df[df["chromhmm_state"] != "UNANNOTATED"].copy()

df["baseMean"] = pd.to_numeric(df["baseMean"])
df["gene_length"] = pd.to_numeric(df["gene_length"])

values = df[["baseMean", "gene_length"]].to_numpy()

if not np.isfinite(values).all():
    raise ValueError("Nonfinite covariates")

if (df["baseMean"] < 0).any():
    raise ValueError("Negative expression")

if (df["gene_length"] <= 0).any():
    raise ValueError("Invalid gene length")

df["log_expression"] = np.log10(df["baseMean"] + 1)
df["log_length"] = np.log10(df["gene_length"])

df["biotype"] = np.select(
    [
        df["gene_type"] == "protein_coding",
        df["gene_type"] == "lncRNA",
    ],
    ["protein_coding", "lncRNA"],
    default="other",
)

background = df[df["signature_group"] == "background"]

# Use the same biotype-specific scaling approach as Phase 14A.2.
scaling = {}

for biotype, pool in background.groupby("biotype"):
    center = pool[FEATURES].mean().to_numpy()
    scale = pool[FEATURES].std(ddof=1).to_numpy()

    if not np.isfinite(scale).all() or (scale <= 0).any():
        raise ValueError(f"Invalid scaling: {biotype}")

    scaling[biotype] = (center, scale)

pairs = []

for group in GROUPS:

    selected = df[df["signature_group"] == group]

    for (chrom, biotype), treated in selected.groupby(
        ["chrom", "biotype"], sort=True
    ):

        controls = background[
            (background["chrom"] == chrom) &
            (background["biotype"] == biotype)
        ]

        if controls.empty:
            continue

        treated = treated.sort_values(
            "gene_id"
        ).reset_index(drop=True)

        controls = controls.sort_values(
            "gene_id"
        ).reset_index(drop=True)

        center, scale = scaling[biotype]

        X = (controls[FEATURES].to_numpy() - center) / scale
        Y = (treated[FEATURES].to_numpy() - center) / scale

        tree = cKDTree(X)

        nearest, _ = tree.query(Y, k=1)
        nearest = np.atleast_1d(nearest)

        # Match the hardest-to-match genes first.
        order = sorted(
            range(len(treated)),
            key=lambda i: (
                -float(nearest[i]),
                str(treated.iloc[i]["gene_id"]),
            ),
        )

        used = set()

        for i in order:

            candidates = tree.query_ball_point(
                Y[i], r=CALIPER
            )

            available = [
                j for j in candidates if j not in used
            ]

            if not available:
                continue

            j = min(
                available,
                key=lambda j: (
                    float(np.linalg.norm(Y[i] - X[j])),
                    str(controls.iloc[j]["gene_id"]),
                ),
            )

            distance = float(np.linalg.norm(Y[i] - X[j]))

            if distance > CALIPER + 1e-10:
                raise RuntimeError("Caliper violation")

            used.add(j)

            a = treated.iloc[i]
            b = controls.iloc[j]

            pairs.append({
                "comparison": group + "_vs_background",
                "chromosome": chrom,
                "biotype": biotype,
                "signature_gene_id": a["gene_id"],
                "background_gene_id": b["gene_id"],
                "signature_state": a["chromhmm_state"],
                "background_state": b["chromhmm_state"],
                "signature_log_expression": a["log_expression"],
                "background_log_expression": b["log_expression"],
                "signature_log_length": a["log_length"],
                "background_log_length": b["log_length"],
                "standardized_distance": distance,
            })

out = pd.DataFrame(pairs)

if out.empty:
    raise RuntimeError("No matches constructed")

if out["standardized_distance"].gt(
    CALIPER + 1e-10
).any():
    raise RuntimeError("Caliper violation")

if out.duplicated(
    ["comparison", "signature_gene_id"]
).any():
    raise RuntimeError("Signature gene duplicated")

if out.duplicated(
    ["comparison", "background_gene_id"]
).any():
    raise RuntimeError("Background gene reused")

gene_chromosome = df.set_index("gene_id")["chrom"]

if not (
    out["signature_gene_id"].map(gene_chromosome).to_numpy()
    == out["chromosome"].to_numpy()
).all():
    raise RuntimeError("Signature chromosome mismatch")

if not (
    out["background_gene_id"].map(gene_chromosome).to_numpy()
    == out["chromosome"].to_numpy()
).all():
    raise RuntimeError("Background chromosome mismatch")

print("=== Phase 14A.5 chromosome-exact matching ===")

for group in GROUPS:

    sub = out[
        out["comparison"] == group + "_vs_background"
    ]

    total = int((df["signature_group"] == group).sum())

    print(f"\n=== {group} ===")
    print("Eligible signature genes:", total)
    print("Matched pairs:", len(sub))
    print("Unmatched:", total - len(sub))
    print(
        "Retention:",
        f"{100 * len(sub) / total:.2f}%"
    )

    if sub.empty:
        raise RuntimeError(f"No matches for {group}")

    print(
        "Maximum matching distance:",
        f"{sub['standardized_distance'].max():.4f}"
    )

    for feature in FEATURES:

        a = sub["signature_" + feature].to_numpy()
        b = sub["background_" + feature].to_numpy()

        pooled_sd = np.sqrt(
            (
                np.var(a, ddof=1) +
                np.var(b, ddof=1)
            ) / 2
        )

        smd = (
            (a.mean() - b.mean()) / pooled_sd
            if pooled_sd > 0 else float("nan")
        )

        print(f"Matched {feature} SMD: {smd:.5f}")

    print(
        "Biotype counts:",
        sub["biotype"].value_counts().to_dict()
    )

    print(
        "Chromosomes represented:",
        sub["chromosome"].nunique()
    )

out = out.sort_values(
    ["comparison", "chromosome", "biotype",
     "signature_gene_id"]
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

out.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    mode="x",
    lineterminator="\n",
)

print("\nPASS: Chromosome-exact matched cohorts constructed")
print("Output:", OUTPUT)
