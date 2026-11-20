#!/usr/bin/env python3
"""Independent, descriptive HepG2 count-direction audit."""

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

OUT = Path("results/phase18")
OUT.mkdir(parents=True, exist_ok=True)

counts_path = Path(
    "results/counts/validation/GSE271000/"
    "GSE271000_validation_gene_counts.tsv"
)

meta_path = Path(
    "metadata/GSE271000_validation_samples.tsv"
)

annotation_path = Path(
    "results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"
)

raw_dir = Path("data/raw/validation/GSE271000")

genes = [
    "SREBF2", "SREBF1", "HMGCR", "HMGCS1",
    "SQLE", "FDFT1", "LDLR", "CYP51A1",
    "DHCR7", "TM7SF2", "MSMO1", "HSD17B7",
    "LSS", "NSDHL", "SC5D", "DHCR24", "EBP"
]

meta = pd.read_csv(meta_path, sep="\t", dtype=str)
df = pd.read_csv(counts_path, sep="\t")

assert len(meta) == 8
assert df.shape == (78932, 9), df.shape
assert df.columns.tolist() == ["gene_id"] + meta.run.tolist()
assert df.gene_id.is_unique

expected = {
    "control": 2,
    "shSREBP2_664": 2,
    "shSREBP2_665": 2,
    "shSREBP2_667": 2
}
assert meta.groupby("group").size().to_dict() == expected

X = df[meta.run.tolist()].to_numpy(dtype=float)

assert np.isfinite(X).all()
assert (X >= 0).all()
assert (X == np.floor(X)).all()

# DESeq2-style median-ratio size-factor estimator using genes
# positive in every sample. This is a descriptive approximation,
# not a fitted negative-binomial DESeq2 contrast.
mask = (X > 0).all(axis=1)
assert mask.sum() > 1000, "Too few fully positive genes"

log_geo = np.mean(np.log(X[mask]), axis=1)
ratios = np.exp(np.log(X[mask]) - log_geo[:, None])
size_factors = np.median(ratios, axis=0)

assert np.isfinite(size_factors).all()
assert (size_factors > 0).all()

normalized = X / size_factors[None, :]

annotations = pd.read_csv(
    annotation_path,
    sep="\t",
    usecols=["gene_id", "gene_symbol"]
)

annotations = annotations[
    annotations.gene_symbol.isin(genes)
].drop_duplicates()

annotations = annotations[
    ~annotations.gene_symbol.duplicated(keep=False)
]

lookup = annotations.set_index("gene_symbol").gene_id

missing = set(genes) - set(lookup.index)
if missing:
    raise ValueError(
        f"Missing or ambiguous gene annotations: {sorted(missing)}"
    )

gene_index = pd.Series(
    np.arange(len(df)), index=df.gene_id
)

records = []

for hp in ["664", "665", "667"]:
    treated = meta.loc[
        meta.group.eq(f"shSREBP2_{hp}"), "run"
    ].tolist()
    controls = meta.loc[
        meta.group.eq("control"), "run"
    ].tolist()

    treated_idx = [
        meta.run.tolist().index(s) for s in treated
    ]
    control_idx = [
        meta.run.tolist().index(s) for s in controls
    ]

    raw = pd.read_csv(
        raw_dir /
        f"GSE271000_DESeq_HepG2_shSREBP2_{hp}.txt.gz",
        sep="\t",
        compression="gzip"
    ).rename(columns={"Unnamed: 0": "gene_symbol"})

    assert raw.gene_symbol.is_unique
    raw = raw.set_index("gene_symbol")

    for gene in genes:
        gid = lookup[gene]
        idx = int(gene_index[gid])

        treatment_mean = normalized[
            idx, treated_idx
        ].mean()
        control_mean = normalized[
            idx, control_idx
        ].mean()

        # Pseudocount only stabilizes near-zero means.
        descriptive_lfc = np.log2(
            (treatment_mean + 0.5) /
            (control_mean + 0.5)
        )

        published_lfc = float(
            raw.loc[gene, "log2FoldChange"]
        )

        records.append({
            "hairpin": hp,
            "gene": gene,
            "gene_id": gid,
            "control_normalized_mean": control_mean,
            "treated_normalized_mean": treatment_mean,
            "count_log2_ratio": descriptive_lfc,
            "deposited_log2FC": published_lfc,
            "same_sign": (
                np.sign(descriptive_lfc) ==
                np.sign(published_lfc)
            ),
            "count_direction_reliable": (
                control_mean >= 10 and treatment_mean >= 10
            )
        })

out = pd.DataFrame(records)
out.to_csv(
    OUT / "phase18e2_count_direction_audit.tsv",
    sep="\t",
    index=False
)

sf = meta[["run", "group"]].copy()
sf["estimated_size_factor"] = size_factors
sf.to_csv(
    OUT / "phase18e2_size_factors.tsv",
    sep="\t",
    index=False
)

print("=== SIZE FACTORS ===")
print(sf.to_string(index=False))

print("\n=== SREBF2 COUNT-LEVEL CHECK ===")
print(
    out[out.gene.eq("SREBF2")][[
        "hairpin", "control_normalized_mean",
        "treated_normalized_mean",
        "count_log2_ratio",
        "deposited_log2FC",
        "same_sign"
    ]].to_string(index=False)
)

print("\n=== DIRECTION CONCORDANCE ===")
summary = out.groupby("hairpin").agg(
    genes=("gene", "size"),
    same_sign=("same_sign", "sum"),
    reliable_count_direction=(
        "count_direction_reliable", "sum"
    )
)
print(summary.to_string())

print("\n=== DISAGREEMENTS ===")
print(
    out[~out.same_sign][[
        "hairpin", "gene",
        "count_log2_ratio", "deposited_log2FC"
    ]].to_string(index=False)
)

if out.same_sign.all():
    print(
        "\nPASS: All selected count-based directions agree "
        "with deposited DESeq2 results."
    )
else:
    print(
        "\nREVIEW: Some directions disagree. "
        "Inspect counts and quantification provenance."
    )

print(
    "\nNOTE: This audit establishes direction consistency "
    "from existing counts; it is not an independent "
    "DESeq2 differential-expression fit."
)
