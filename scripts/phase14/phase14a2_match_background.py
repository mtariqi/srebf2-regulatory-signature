#!/usr/bin/env python3
"""Phase 14A.2: deterministic biotype-exact background matching.

Greedy nearest-neighbor matching without replacement.
Diagnostic matching design, not globally optimal assignment.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SOURCE = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)
OUTPUT = Path(
    "results/phase14/SREBF2_PHASE14A2_MATCHED_PAIRS.tsv"
)
CALIPER = 0.25
FEATURES = ["log_expression", "log_length"]
GROUPS = ["increased_in_KO", "decreased_in_KO"]
BIOTYPES = ["protein_coding", "lncRNA", "other"]

if OUTPUT.exists():
    raise FileExistsError(f"Refusing overwrite: {OUTPUT}")

df = pd.read_csv(SOURCE, sep="\t")
assert df["gene_id"].is_unique
assert len(df) == 17097

df = df[df["chromhmm_state"] != "UNANNOTATED"].copy()
df["baseMean"] = pd.to_numeric(df["baseMean"])
df["gene_length"] = pd.to_numeric(df["gene_length"])

if not np.isfinite(df[["baseMean", "gene_length"]]).all().all():
    raise ValueError("Nonfinite covariates")
if (df["baseMean"] < 0).any() or (df["gene_length"] <= 0).any():
    raise ValueError("Invalid covariates")

df["log_expression"] = np.log10(df["baseMean"] + 1)
df["log_length"] = np.log10(df["gene_length"])
df["biotype"] = np.select(
    [df["gene_type"].eq("protein_coding"),
     df["gene_type"].eq("lncRNA")],
    ["protein_coding", "lncRNA"],
    default="other",
)

background = df[df["signature_group"] == "background"]
pairs = []

for group in GROUPS:
    for biotype in BIOTYPES:
        treated = (
            df[(df["signature_group"] == group) &
               (df["biotype"] == biotype)]
            .sort_values("gene_id")
            .reset_index(drop=True)
        )
        controls = (
            background[background["biotype"] == biotype]
            .sort_values("gene_id")
            .reset_index(drop=True)
        )

        center = controls[FEATURES].mean().to_numpy()
        scale = controls[FEATURES].std(ddof=1).to_numpy()

        if np.any(scale <= 0) or not np.isfinite(scale).all():
            raise ValueError(f"Invalid scaling: {biotype}")

        X = (controls[FEATURES].to_numpy() - center) / scale
        Y = (treated[FEATURES].to_numpy() - center) / scale
        tree = cKDTree(X)

        nearest, _ = tree.query(Y, k=1)

        # Handle the hardest-to-match genes first.
        order = sorted(
            range(len(treated)),
            key=lambda i: (-float(nearest[i]),
                           treated.iloc[i]["gene_id"])
        )

        used = set()

        for i in order:
            available = tree.query_ball_point(
                Y[i], r=CALIPER
            )

            available = [j for j in available if j not in used]

            if not available:
                continue

            j = min(
                available,
                key=lambda j: (
                    float(np.linalg.norm(Y[i] - X[j])),
                    controls.iloc[j]["gene_id"]
                ),
            )

            used.add(j)

            a = treated.iloc[i]
            b = controls.iloc[j]
            distance = float(np.linalg.norm(Y[i] - X[j]))

            pairs.append({
                "comparison": group + "_vs_background",
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
    raise RuntimeError("No valid matches found")

if out["standardized_distance"].gt(CALIPER + 1e-10).any():
    raise RuntimeError("Matching caliper violation")

for group in GROUPS:
    part = out[out["comparison"] == group + "_vs_background"]

    if part["signature_gene_id"].duplicated().any():
        raise RuntimeError("Signature gene matched twice")

    if part["background_gene_id"].duplicated().any():
        raise RuntimeError("Background gene reused in a comparison")

    original_n = int((df["signature_group"] == group).sum())

    print(f"\n=== {group} ===")
    print(f"Signature genes: {original_n}")
    print(f"Matched pairs: {len(part)}")
    print(f"Unmatched genes: {original_n - len(part)}")
    print(
        "Retention:",
        f"{100 * len(part) / original_n:.2f}%"
    )
    print(
        "Maximum distance:",
        f"{part['standardized_distance'].max():.4f}"
    )

    for feature in FEATURES:
        a = part["signature_" + feature].to_numpy()
        b = part["background_" + feature].to_numpy()

        pooled_sd = np.sqrt(
            (np.var(a, ddof=1) + np.var(b, ddof=1)) / 2
        )

        smd = ((a.mean() - b.mean()) / pooled_sd
               if pooled_sd > 0 else np.nan)

        print(f"Matched {feature} SMD: {smd:.5f}")

    print(
        "Matched biotypes:",
        part["biotype"].value_counts().to_dict()
    )

out = out.sort_values(
    ["comparison", "biotype", "signature_gene_id"]
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    mode="x",
    lineterminator="\n",
)

print("\nPASS: Matched cohorts constructed")
print("Output:", OUTPUT)
