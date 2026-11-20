# Phase 19 — SREBF2 statistical reproducibility audit

**Status:** numerical checks completed from archived files and user-reported Ubuntu audit output on 2026-11-20. Regional block tables have been produced locally and require synchronization to GitHub.

## Study context

This audit concerns the integrative reanalysis of HeLa SREBF2-knockout RNA-seq and historical Roadmap E117 HeLa-S3 chromatin-state annotations. The chromatin analysis is an exploratory *reference-state association* analysis, **not** direct binding, chromatin remodeling, or experimental replication.

## A. Initial E117 state-screening audit

- 45 Fisher exact tests (15 ChromHMM states × 3 contrasts)
- Source: `results/phase13/SREBF2_E117_STATE_FISHER_ENRICHMENT.tsv`
- Independently recalculated Fisher exact p-values: all matched
- Independently recalculated odds ratios: all matched
- Recomputed BH FDR against stored `fdr_bh_45tests`: maximum absolute difference **0.0**
- Significant at BH FDR < 0.05: **12 tests**
- Increased-in-KO versus background, historical `15_Quies`: unmatched OR **3.244201**, BH FDR **4.169555e-28**
- Increased-in-KO versus background, `1_TssA`: unmatched OR **0.357092**, BH FDR **1.005741e-28**

The subsequent six selected chromatin-state follow-ups are **post-screen**, not independent prespecified confirmatory tests; BH adjustment within each follow-up family does not account for the entire selection procedure.

## B. Follow-up adjusted models

- Source: `results/phase14/SREBF2_PHASE14A1_MODEL_AUDIT.tsv`
- 6 selected states × 3 contrasts = **18 models**
- **18/18 converged**; all full-rank; all numerical estimates finite
- Refit versus frozen odds-ratio maximum absolute difference: **4.440892098500626e-16**
- Maximum relative difference: **2.1543498548704486e-16**

## C. Chromosome-exact matching

- Source: `results/phase14/SREBF2_PHASE14A5_CHROMOSOME_EXACT_PAIRS.tsv`
- Published analysis script: `scripts/phase14/phase14a6_chromosome_exact_enrichment.py`
- Primary increased-in-KO × `15_Quies`: **425 pairs**, 25 both-state, **106** signature-only, **55** background-only, 239 neither; matched OR **1.9272727273**, exact 95% CI **1.3787896781–2.7202333209**, exact two-sided p **7.149120725847529e-05**, BH across 2 matched tests **0.00014298241451695057**
- Secondary decreased-in-KO × `14_ReprPCWk`: **677 pairs**, 0 both-state, 25 signature-only, 20 background-only, 632 neither; matched OR **1.25**, exact 95% CI **0.6665930488–2.3733192671**, exact p **0.5514843298025198**
- Both comparisons were independently reconstructed from pair-level records.

## D. Chromosomal sensitivity

- Chromosome-level net-discordance signs in the 425-pair quiescent comparison: **16 positive, 2 negative, 5 zero** (23 chromosomes)
- Exhaustive two-sided sign-flip enumeration for the 18 nonzero chromosome contributions: **546 / 262,144 = 0.00208282470703125**, assuming chromosome-level sign symmetry (alternative to earlier Monte Carlo estimate).
- Local chromosome-exact leave-one-chromosome-out run: **23 exclusions**, all OR > 1, OR range **1.781818–2.102041**, maximum nominal exact p **0.000833**
- Local table: `results/phase14/SREBF2_PHASE19_CHROMOSOME_EXACT_LOCO_AUDIT.tsv` (verify this file has been pushed).

## E. 1-Mb regional-block sensitivity

Using 0-based signature TSS coordinates (`tss_0based`) and chromosome from `results/phase13/SREBF2_E117_TSS_CHROMATIN_STATES.tsv`:

- Primary chromosome-exact quiescent cohort: **425 matched pairs**
- **340** unique 1-Mb signature-TSS blocks; **141** informative blocks
- Overall discordance: 106 signature-only vs 55 background-only (net +51)
- Largest individual net discordance contribution: **+3**, attained by `chr12:9`, `chr19:23`, `chr19:21`
- Leave-one-block-out OR range **1.8727272727–1.9811320755**
- Largest exact two-sided p: **0.00016616729130786213**
- All leave-one-block-out OR > 1; all nominal p < 0.05
- Locally produced tables **not yet confirmed uploaded**:
  - `results/phase19/SREBF2_PHASE19_1MB_BLOCK_SUMMARY.tsv`
  - `results/phase19/SREBF2_PHASE19_LEAVE_ONE_BLOCK_OUT.tsv`

LOCO and block-exclusion estimates are overlapping sensitivity analyses; their nominal p-values should not be interpreted as independent confirmations. Block-level inference may still depend on spatial correlation and matching design.

## F. GSEA numerical sensitivity

The primary and `fgseaMultilevel(eps=0)` results were compared on the same named cholesterol/sterol pathways:

| Pathway | Original NES | eps=0 NES | Original BH FDR | eps=0 BH FDR |
|---|---:|---:|---:|---:|
| GO sterol biosynthetic process | -2.9484996375 | -2.9740972129 | 4.046153846e-08 | 4.517380583e-15 |
| GO cholesterol biosynthetic process | -2.9100705622 | -2.9188418442 | 4.046153846e-08 | 1.615209881e-13 |
| Reactome cholesterol biosynthesis | -3.0277531432 | -3.0483608945 | 6.060869565e-09 | 5.036248071e-19 |

Source tables contained 5,260 GO and 1,394 Reactome pathway rows, with no malformed rows. `Figure_GSEA_Running.png` and `Figure_GSEA_Running_eps0.png` have identical SHA-256 checksums `8857023def37432fc4c2627c8bad4f512f17556119aa4fa15d495429ae721251`; these are duplicate *visualizations*, not independent figures. Preserve original figure labels; report eps=0 differences numerically in a supplementary table.

## Remaining reporting tasks

1. Push local Phase 19 TSVs from Ubuntu after inspecting `git status`.
2. Attach a reproducible single-command audit script / environment lock and hash manifests.
3. Review the documented timing and rationale for selecting six follow-up states after the 45-test screen.
4. Consider dependence-aware block/bootstrap or chromosome-cluster inference as an additional sensitivity; do not describe existing tests as independent validation.
5. Update the final manuscript and supplement, including original versus eps=0 statistics in separate columns.
