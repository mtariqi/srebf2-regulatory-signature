# SREBF2 PCA and exploratory GNN extension

**Scientific purpose:** Characterize experimental sample structure in HeLa RNA-seq and explore the topology of a broader STRING functional-association graph. These are **separate analyses**: RNA-seq PCA visualizes experimental samples; GNN embedding plots visualize proteins. Neither demonstrates causality, direct SREBP2 binding, or chromatin-state transitions.

## Verified inputs

- `results/counts/GSE267018_raw_count_matrix.tsv` — original **12-sample** HeLa RNA-seq counts, already tracked in the SREBF2 repository.
- `metadata/GSE267018_analysis_samples.tsv` — 12 runs, 2 genotypes × 2 treatments, 3 biological replicates per group.
- `srebf2_pca_gnn/data/hela_sterol_core_effects.tsv` — 12 original HeLa knockout log2FC estimates, **used only as visual overlays**, not model labels.
- `srebf2_pca_gnn/data/string_core_5_verified.tsv` — 10 high-scoring STRING functional associations between five previously selected sterol enzymes. **The five-node complete graph cannot be used to assess GNN predictive ability.**

## Installation (Ubuntu; project venv)

```bash
cd ~/srebf2-regulatory-signature
source .venv-srebf2-validation/bin/activate
python -m pip install scikit-learn torch requests pandas numpy networkx matplotlib
```

To use recommended DESeq2 VST PCA, separately ensure R libraries `DESeq2` and `matrixStats` are available (do not install system-wide packages casually).

## 1. Run sample-level PCA using original HeLa counts

```bash
python srebf2_pca_gnn/scripts/01_pca_samples.py \
  --counts results/counts/GSE267018_raw_count_matrix.tsv \
  --metadata metadata/GSE267018_analysis_samples.tsv \
  --out srebf2_pca_gnn/figures
```

This generates **both** `PCA_all12.png` and `PCA_exclude_SRR28966285.png` with score and explained-variance TSVs. Python PCA uses library-size normalized log2(CPM+1), highest-variance expressed genes, and PCA on **samples**, not genes. This is exploratory QC; the alternative **preferred** DESeq2 variance-stabilized approach is:

```bash
Rscript srebf2_pca_gnn/scripts/01b_pca_deseq2_vst.R \
  results/counts/GSE267018_raw_count_matrix.tsv \
  metadata/GSE267018_analysis_samples.tsv \
  srebf2_pca_gnn/figures
```

Report transformation and filtering parameters transparently. Excluding run SRR28966285 is a sensitivity analysis, **not** an automatic basis for removing it from the primary model.

## 2. Expand the STRING human network

```bash
python srebf2_pca_gnn/scripts/02_fetch_string_expanded.py \
  --out srebf2_pca_gnn/data \
  --add-nodes 70 \
  --score 700
```

This requests the 12 core enzymes and 70 additional neighbors using the public STRING API. It records the API request, response, network type (functional), score threshold, taxon 9606, and retrieval time. Available expansion settings and behavior may differ by STRING release; review returned sizes and edge coverage.

## 3. Train and evaluate exploratory GNN

```bash
python srebf2_pca_gnn/scripts/03_gnn_string_autoencoder.py \
  --edges srebf2_pca_gnn/data/string_expanded_edges.tsv \
  --effects srebf2_pca_gnn/data/hela_sterol_core_effects.tsv \
  --out srebf2_pca_gnn/figures
```

The two-layer GCN graph autoencoder learns 16-dimensional protein embeddings and predicts whether **withheld STRING association edges** are present. Edge splits preserve a spanning forest for connectivity; negatives are nonedges in the full retrieved graph; structural features are computed **only from the training graph**, preventing trivial test-edge leakage. Evaluation includes ROC AUC and average precision against a common-neighbors heuristic on the same split. Test performance measures reconstruction of STRING associations, **not** SREBF2 knockout response prediction or biological discovery. RNA-seq fold changes are visual overlays only, never supervision.

**Fail-safe:** The script refuses to train on a graph with fewer than 30 nodes, a complete graph, or insufficient independent held-out positive/negative edges. The verified five-gene network by itself does not satisfy these requirements. Do not interpret a single held-out split as proof of GNN generalizability.

## Outputs and reporting policy

| Figure | Source | Scientific interpretation |
|---|---|---|
| RNA PCA, 12 samples | GSE267018 count matrix; sample metadata | Exploratory sample relationships; not differential-expression significance |
| PCA leave-one-out sensitivity | Same counts excluding SRR28966285 | Robustness of descriptive sample geometry |
| GNN embedding view | Expanded STRING associations | Exploratory functional-association topology; HeLa log2FC visual overlay |

- The PCA outputs **must be generated locally from actual counts**; none are synthesized for this package.
- The GNN figure **must be generated locally from an expanded network**; the 10-edge complete five-gene network is insufficient.
- STRING combines evidence from databases, coexpression, text mining, and experiments. Functional associations are not confirmed direct physical interactions.
- Protein-neighbor sampling, pathway-based seed selection, graph topology, and database annotation introduce substantial selection bias. Nonedges are unlabeled, **not confirmed negative interactions**.
- Test several network thresholds, graph expansions, and randomized splits before claiming robustness. Comparing to trivial baselines is required; performance alone does not demonstrate regulatory mechanisms.
- A stronger **future** model for KO responsiveness would require independent perturbation labels, non-leaking features, gene-level train/test separation, and external biological validation.

## References

- [DESeq2 guide: normalization, VST and PCA](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [STRING API documentation](https://string-db.org/help/api/)
