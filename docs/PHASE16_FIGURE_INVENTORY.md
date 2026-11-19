# Phase 16 figure inventory and provenance

This repository reanalyzes public HeLa **SREBF2** knockout RNA-seq and uses historical reference chromatin annotations. The Phase 16 figures are graphical summaries of the frozen analyses, **not independent experiments**.

## Latest figure set (local output paths)

| Figure | Description | Script / evidence |
|---|---|---|
| `Figure_Volcano_MA.png` and `.pdf` | KO-versus-WT volcano and MA plots, 12 sterol-core labels | `scripts/phase16/phase16a_volcano_ma.R` (generated in user's Ubuntu checkout) |
| `Figure_GSEA_Running.png` and `.pdf` | Genuine GO sterol, GO cholesterol, Reactome cholesterol GSEA running curves | `scripts/phase16/phase16b_gsea_running_plots.R` |
| `Figure_GSEA_Running_eps0.png` and `.pdf` | Numerical `fgseaMultilevel(eps=0)` sensitivity; preserve original GSEA statistics separately | Same script with `--eps0` |
| `Figure_SterolCore_Heatmap.png` and `.pdf` | Twelve-gene VST row-Z-score heatmap of 12 samples | `scripts/phase16/phase16c_sterol_heatmap.R` |
| `Figure_Metabolic_Programs.png` and `.pdf` | Negative NES comparison and *descriptive* respiratory-complex schematic | `scripts/phase16/phase16e_metabolic_programs.py`; `results/phase16/figures/Figure_Metabolic_Programs_source.tsv` |

**Important upload status:** This inventory documents figure files available in the local analysis or manuscript uploads. Listing a file here does **not** mean its binary image has been committed to GitHub. Verify file presence in the repository independently.

### Reproduce the metabolic-programs figure

On Ubuntu, use a virtual environment (PEP 668-compatible):

```bash
cd ~/srebf2-regulatory-signature
python3 -m venv .venv-phase16
source .venv-phase16/bin/activate
python -m pip install matplotlib numpy
python scripts/phase16/phase16e_metabolic_programs.py \
  --table results/phase16/figures/Figure_Metabolic_Programs_source.tsv \
  --output-prefix results/phase16/figures/Figure_Metabolic_Programs
```

### Publish all latest graphics from the analysis workstation

```bash
cd ~/srebf2-regulatory-signature
git pull --ff-only origin main
git add scripts/phase16/ results/phase16/figures/
git diff --cached --stat
git commit -m "Publish Phase 16 original and metabolic-programs figures"
git push origin main
```

Before staging, review the files and any sensitive/raw material. Do not use `git add .`. If a file is ignored, check the ignore rule first. If the working branch differs or `git pull` reports conflicts, resolve those before proceeding.

### Limitations

The negative enrichment statistics for metabolic programs are frozen *manuscript-reported* values, not independently recomputed by the Python visualization. GO and Reactome sets overlap; findings are not biological replication. The 57 oxidative-phosphorylation leading-edge genes and respiratory-chain families are described at group level; the diagram is **not** an empirical per-complex measurement. No oxygen consumption, respiratory activity, ATP production, sterol flux or direct chromatin transition was measured by this figure. Exact Reactome translation values require a separate verified export.
