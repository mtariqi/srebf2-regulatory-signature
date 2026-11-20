#!/usr/bin/env Rscript

# Phase 17C: contrast-level DoRothEA/VIPER analysis
# Signed DESeq2 Wald statistics: KO versus WT, depletion
# Outputs are exploratory, not independent validation.

phase17_lib <- path.expand("~/R/srebf2-phase17-library")
if (!dir.exists(phase17_lib)) {
  stop("Phase 17 R library is missing")
}
.libPaths(c(phase17_lib, .libPaths()))

suppressPackageStartupMessages({
  library(dorothea)
  library(decoupleR)
})

outdir <- "results/phase17"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ranked <- read.delim(
  file.path(outdir, "phase17b_ranked_unique_symbols.tsv"),
  check.names = FALSE
)

stopifnot(
  all(c("gene_symbol", "stat") %in% colnames(ranked)),
  !anyDuplicated(ranked$gene_symbol),
  all(is.finite(ranked$stat))
)

# Matrix: genes × one signed DE contrast.
mat <- matrix(
  ranked$stat,
  ncol = 1,
  dimnames = list(
    ranked$gene_symbol,
    "KO_vs_WT_depletion"
  )
)

data("dorothea_hs", package = "dorothea")

regulons <- as.data.frame(dorothea_hs)

stopifnot(
  all(c("tf", "target", "mor", "confidence") %in%
      colnames(regulons))
)

# Primary analysis: evidence confidence A-C.
regulons <- regulons[
  regulons$confidence %in% c("A", "B", "C") &
  regulons$target %in% rownames(mat),
]

network <- data.frame(
  source = regulons$tf,
  target = regulons$target,
  mor = as.numeric(regulons$mor)
)

network <- network[
  is.finite(network$mor) &
  !is.na(network$source) &
  !is.na(network$target),
]

# TF-target pairs should be unique in the reference.
if (anyDuplicated(network[c("source", "target")])) {
  stop("Duplicate TF-target edges detected; inspect reference")
}

cat("=== PHASE 17C: VIPER PRIMARY ANALYSIS ===\n")
cat("Ranked genes:", nrow(mat), "\n")
cat("Regulon edges:", nrow(network), "\n")
cat("Candidate TFs:", length(unique(network$source)), "\n")

# Preserve input Wald statistics using method="none".
# minsize = 10 ensures adequate target overlap.
results <- decoupleR::run_viper(
  mat = mat,
  network = network,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 10,
  method = "none",
  pleiotropy = FALSE,
  verbose = FALSE
)

results <- as.data.frame(results)

stopifnot(
  nrow(results) > 0,
  all(c("source", "score") %in% colnames(results)),
  !anyDuplicated(results$source),
  all(is.finite(results$score))
)

# Report modeled target coverage.
coverage <- read.delim(
  file.path(outdir, "phase17b_regulon_coverage.tsv")
)

results <- merge(
  results,
  coverage[, c(
    "tf", "total_targets",
    "measured_targets", "coverage_fraction"
  )],
  by.x = "source",
  by.y = "tf",
  all.x = TRUE
)

# Nominal p-values, if returned, are model-derived
# rather than independent sample-level validation.
if ("p_value" %in% names(results)) {
  results$fdr_bh <- p.adjust(
    results$p_value,
    method = "BH"
  )
}

results <- results[
  order(-abs(results$score)),
]

write.table(
  results,
  file.path(outdir, "phase17c_viper_activity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n=== TF ACTIVITY SUMMARY ===\n")
cat("TFs scored:", nrow(results), "\n")

cat("\n=== CANDIDATE REGULATORS ===\n")
candidates <- c("SREBF2", "SREBF1", "MYC", "SP1")

candidate_results <- results[
  results$source %in% candidates,
]

print(
  candidate_results[
    order(match(candidate_results$source, candidates)),
  ],
  row.names = FALSE
)

cat("\n=== TOP 20 TFs BY ABSOLUTE ACTIVITY ===\n")

print(
  head(results, 20),
  row.names = FALSE
)

cat("\n=== PACKAGE VERSIONS ===\n")
cat("dorothea:", as.character(packageVersion("dorothea")), "\n")
cat("decoupleR:", as.character(packageVersion("decoupleR")), "\n")
cat("viper:", as.character(packageVersion("viper")), "\n")

cat("\nPASS: Phase 17C VIPER analysis complete\n")
