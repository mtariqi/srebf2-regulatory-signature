#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p results/phase16/figures results/phase16/tables
# The actual all12 DESeq2 file is tracked in the repository.
DE='results/de/discovery_v1/all12/KO_vs_WT_depletion.tsv'
COUNT='results/counts/GSE267018_raw_count_matrix.tsv'
META='metadata/GSE267018_analysis_samples.tsv'
GO='results/phase15/phase15e_enrichment/GO_BP_GSEA.rds'
REACT='results/phase15/phase15e_enrichment/Reactome_GSEA.rds'
for file in "$DE" "$COUNT" "$META" "$GO" "$REACT"; do
  if [[ ! -f "$file" ]]; then echo "ERROR: Missing $file" >&2; exit 2; fi
done
Rscript scripts/phase16/phase16a_volcano_ma.R "$DE" results/phase16/figures/Figure_Volcano_MA
Rscript scripts/phase16/phase16b_gsea_running_plots.R "$GO" "$REACT" results/phase16/figures/Figure_GSEA_Running
Rscript scripts/phase16/phase16c_sterol_heatmap.R "$COUNT" "$META" results/phase16/figures/Figure_SterolCore_Heatmap
sha256sum "$DE" "$COUNT" "$META" "$GO" "$REACT" > results/phase16/tables/input_SHA256.txt
Rscript -e 'writeLines(capture.output(sessionInfo()), "results/phase16/tables/R_sessionInfo.txt")'
echo 'PASS: Three Phase 16 figure sets generated from actual repository inputs.'
