#!/usr/bin/env python3
"""Phase 14A.1: numerical audit of frozen Phase 13H models."""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

SOURCE = Path(
    "results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv"
)

REFERENCE = Path(
    "results/phase13/SREBF2_E117_SPLINE_SENSITIVITY_v2.tsv"
)

OUTPUT = Path(
    "results/phase14/SREBF2_PHASE14A1_MODEL_AUDIT.tsv"
)

STATES = [
    "1_TssA",
    "2_TssAFlnk",
    "5_TxWk",
    "7_Enh",
    "14_ReprPCWk",
    "15_Quies",
]

CONTRASTS = [
    ("decreased_in_KO", "background"),
    ("increased_in_KO", "background"),
    ("increased_in_KO", "decreased_in_KO"),
]

FORMULA = (
    'outcome ~ state_present '
    '+ cr(log_baseMean, df=4, constraints="center") '
    '+ cr(log_gene_length, df=4, constraints="center") '
    '+ C(biotype_group)'
)

if OUTPUT.exists():
    raise FileExistsError(
        f"Refusing to overwrite existing file: {OUTPUT}"
    )

df = pd.read_csv(SOURCE, sep="\t")
frozen = pd.read_csv(REFERENCE, sep="\t")

if df["gene_id"].duplicated().any():
    raise ValueError("Duplicate gene IDs")

if len(df) != 17097:
    raise ValueError("Unexpected input gene count")

if len(frozen) != 18:
    raise ValueError("Expected 18 frozen model results")

df = df.loc[
    df["chromhmm_state"] != "UNANNOTATED"
].copy()

df["baseMean"] = pd.to_numeric(
    df["baseMean"], errors="raise"
)

df["gene_length"] = pd.to_numeric(
    df["gene_length"], errors="raise"
)

if not np.isfinite(
    df[["baseMean", "gene_length"]].to_numpy()
).all():
    raise ValueError("Nonfinite covariates")

if (df["baseMean"] < 0).any():
    raise ValueError("Negative expression")

if (df["gene_length"] <= 0).any():
    raise ValueError("Invalid gene lengths")

df["log_baseMean"] = np.log10(df["baseMean"] + 1)
df["log_gene_length"] = np.log10(df["gene_length"])

df["biotype_group"] = np.select(
    [
        df["gene_type"] == "protein_coding",
        df["gene_type"] == "lncRNA",
    ],
    ["protein_coding", "lncRNA"],
    default="other",
)

results = []

for group_a, group_b in CONTRASTS:

    sub = df.loc[
        df["signature_group"].isin([group_a, group_b])
    ].copy()

    sub["outcome"] = (
        sub["signature_group"] == group_a
    ).astype(int)

    for state in STATES:

        data = sub.copy()

        data["state_present"] = (
            data["chromhmm_state"] == state
        ).astype(int)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")

            fit = smf.glm(
                FORMULA,
                data=data,
                family=sm.families.Binomial(),
            ).fit(maxiter=200)

        X = fit.model.exog

        rank = int(np.linalg.matrix_rank(X))
        n_columns = int(X.shape[1])
        condition = float(np.linalg.cond(X))

        beta = float(fit.params["state_present"])
        se = float(fit.bse["state_present"])
        calculated_or = float(np.exp(beta))

        contrast = f"{group_a}_vs_{group_b}"

        match = frozen.loc[
            (frozen["contrast"] == contrast) &
            (frozen["state"] == state)
        ]

        if len(match) != 1:
            raise ValueError(
                f"Missing or duplicate frozen result: "
                f"{contrast}, {state}"
            )

        frozen_or = float(match.iloc[0]["adjusted_or"])

        absolute_difference = abs(
            calculated_or - frozen_or
        )

        relative_difference = (
            absolute_difference / frozen_or
        )

        warning_text = "; ".join(
            sorted(set(str(w.message) for w in caught))
        )

        numeric_ok = bool(
            np.isfinite([
                condition,
                beta,
                se,
                calculated_or,
                float(fit.pvalues["state_present"]),
            ]).all()
        )

        results.append({
            "contrast": contrast,
            "state": state,
            "n_observations": int(fit.nobs),
            "converged": bool(fit.converged),
            "iterations": int(
                fit.fit_history["iteration"]
            ),
            "matrix_columns": n_columns,
            "matrix_rank": rank,
            "full_column_rank": rank == n_columns,
            "condition_number": condition,
            "beta": beta,
            "standard_error": se,
            "calculated_or": calculated_or,
            "frozen_or": frozen_or,
            "absolute_or_difference": absolute_difference,
            "relative_or_difference": relative_difference,
            "numeric_finite": numeric_ok,
            "warnings": warning_text,
        })

out = pd.DataFrame(results)

if len(out) != 18:
    raise RuntimeError("Expected 18 models")

if out.duplicated(["contrast", "state"]).any():
    raise RuntimeError("Duplicate audit records")

# Require close agreement with the frozen Phase 13H
# odds ratios, rather than accepting rounded comparisons.
reproduced = bool(
    np.allclose(
        out["calculated_or"],
        out["frozen_or"],
        rtol=1e-8,
        atol=1e-10,
    )
)

out.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    mode="x",
    lineterminator="\n",
)

print("=== Phase 14A.1 numerical stability audit ===")
print("Models:", len(out))
print("Converged:", int(out["converged"].sum()))
print(
    "Full-rank models:",
    int(out["full_column_rank"].sum())
)
print(
    "Numerically finite models:",
    int(out["numeric_finite"].sum())
)
print(
    "Models with warnings:",
    int(out["warnings"].fillna("").str.len().gt(0).sum())
)
print(
    "Maximum condition number:",
    float(out["condition_number"].max())
)
print(
    "Maximum OR difference:",
    float(out["absolute_or_difference"].max())
)
print(
    "Frozen odds ratios reproduced:",
    reproduced
)

print("\n=== Per-model diagnostics ===")
print(
    out[
        [
            "contrast",
            "state",
            "converged",
            "full_column_rank",
            "condition_number",
            "relative_or_difference",
        ]
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.5g}",
    )
)

print("\nOutput:", OUTPUT)

if not (
    out["converged"].all()
    and out["full_column_rank"].all()
    and out["numeric_finite"].all()
    and out["warnings"].fillna("").eq("").all()
    and reproduced
):
    raise RuntimeError(
        "Phase 14A.1 audit failed; inspect the output"
    )

print("\nPASS: Phase 14A.1 consistency and rank checks")
