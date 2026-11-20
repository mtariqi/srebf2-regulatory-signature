#!/usr/bin/env python3
"""Phase 18F: Cross-context sterol-core validation figure."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path("results/phase18")
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

source = pd.read_csv(
    ROOT / "phase18c_sterol_core_gene_level.tsv",
    sep="\t"
)

audit = pd.read_csv(
    ROOT / "phase18e2_count_direction_audit.tsv",
    sep="\t"
)

genes = [
    "TM7SF2", "CYP51A1", "MSMO1", "SQLE",
    "FDFT1", "HSD17B7", "LSS", "NSDHL",
    "DHCR7", "EBP", "SC5D", "DHCR24"
]

comparisons = [
    "HeLa_discovery",
    "HAP1_minus_FBS",
    "HAP1_plus_FBS",
    "HepG2_shRNA_664",
    "HepG2_shRNA_665",
    "HepG2_shRNA_667"
]

# Ensure the source contains exactly the expected core.
assert set(source.gene_symbol) == set(genes)

check = source.groupby(
    ["gene_symbol", "comparison"]
).size()

assert (check == 1).all()
assert len(check) == 60

# Confirm HepG2 count-level signs agree with deposited signs.
sub = audit[audit.gene.isin(genes)]
assert len(sub) == 36
assert sub.same_sign.fillna(False).all()

discovery = (
    source.drop_duplicates("gene_symbol")
    [["gene_symbol", "discovery_log2FC", "discovery_padj"]]
    .rename(columns={
        "discovery_log2FC": "log2FC",
        "discovery_padj": "padj"
    })
)
discovery["comparison"] = "HeLa_discovery"

validation = source[[
    "gene_symbol", "comparison",
    "validation_log2FC", "validation_padj"
]].rename(columns={
    "validation_log2FC": "log2FC",
    "validation_padj": "padj"
})

combined = pd.concat(
    [discovery, validation], ignore_index=True
)

fc = combined.pivot(
    index="gene_symbol",
    columns="comparison",
    values="log2FC"
).reindex(index=genes, columns=comparisons)

fdr = combined.pivot(
    index="gene_symbol",
    columns="comparison",
    values="padj"
).reindex(index=genes, columns=comparisons)

assert fc.shape == (12, 6)
assert fc.notna().all().all()
assert (fc < 0).all().all(), (
    "Some fold changes are not negative; "
    "inspect rather than suppress them."
)

assert fdr.shape == fc.shape

# Keep the full source values in the accompanying TSV.
combined.to_csv(
    ROOT / "phase18f_all_contrast_gene_values.tsv",
    sep="\t",
    index=False
)

fig, ax = plt.subplots(figsize=(12.5, 8))

# Symmetric scale preserves the meaning of zero.
limit = max(3.5, float(np.abs(fc.values).max()))
norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

img = ax.imshow(
    fc.values,
    cmap="RdBu_r",
    norm=norm,
    aspect="auto"
)

labels = [
    "HeLa\nKO / WT",
    "HAP1\n−FBS",
    "HAP1\n+FBS",
    "HepG2\nsh664",
    "HepG2\nsh665",
    "HepG2\nsh667"
]

ax.set_xticks(np.arange(6))
ax.set_xticklabels(labels, fontsize=10)
ax.set_yticks(np.arange(12))
ax.set_yticklabels(genes, fontsize=11)

ax.tick_params(length=0)
ax.set_xticks(np.arange(-.5, 6, 1), minor=True)
ax.set_yticks(np.arange(-.5, 12, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.5)
ax.tick_params(which="minor", bottom=False, left=False)

for i, gene in enumerate(genes):
    for j, comp in enumerate(comparisons):
        value = fc.loc[gene, comp]
        padj = fdr.loc[gene, comp]

        star = "*" if pd.notna(padj) and padj < .05 else ""

        color = (
            "white" if abs(value) > limit * .55
            else "#222222"
        )

        ax.text(
            j, i, f"{value:.2f}{star}",
            ha="center",
            va="center",
            fontsize=9,
            color=color,
            fontweight="normal"
        )

ax.axvline(0.5, color="#333333", linewidth=2)
ax.axvline(2.5, color="#333333", linewidth=2)

cbar = fig.colorbar(
    img, ax=ax, fraction=.035, pad=.03
)
cbar.set_label("Log2 fold change", fontsize=10)

ax.set_title(
    "Cross-context SREBF2 perturbation:\n"
    "12-gene sterol-biosynthesis signature",
    fontsize=15, fontweight="bold", pad=20
)

fig.text(
    .12, .025,
    "* BH-adjusted p < 0.05 within the respective "
    "source contrast.\n"
    "HepG2 direction corroborated using existing "
    "normalized gene counts; source DESeq2 contrast "
    "command not recovered.",
    fontsize=9
)

fig.subplots_adjust(
    left=.16, right=.90, top=.87, bottom=.16
)

for ext in ["png", "pdf"]:
    fig.savefig(
        OUT / f"Figure_Phase18_Sterol_Core_Validation.{ext}",
        dpi=300,
        bbox_inches="tight"
    )

plt.close(fig)

print("=== PHASE 18F RESULTS ===")
print("Genes:", fc.shape[0])
print("Comparisons:", fc.shape[1])
print("Negative fold changes:", int((fc < 0).sum().sum()))
print("FDR < 0.05:", int((fdr < .05).sum().sum()))
print("HepG2 count-level direction audit: PASS")
print("Saved PNG, PDF, and source TSV.")
print("PASS: Phase 18F completed")
