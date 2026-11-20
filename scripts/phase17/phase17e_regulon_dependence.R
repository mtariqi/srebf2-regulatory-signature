#!/usr/bin/env Rscript

lib <- path.expand("~/R/srebf2-phase17-library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(dorothea)
  library(decoupleR)
})

outdir <- "results/phase17"

ranked <- read.delim(
  file.path(outdir, "phase17b_ranked_unique_symbols.tsv")
)

stopifnot(
  !anyDuplicated(ranked$gene_symbol),
  all(is.finite(ranked$stat))
)

data("dorothea_hs", package = "dorothea")

reg <- as.data.frame(dorothea_hs)
reg <- reg[
  reg$confidence %in% c("A", "B", "C"),
]

network <- data.frame(
  source = reg$tf,
  target = reg$target,
  mor = as.numeric(reg$mor)
)

network <- network[
  is.finite(network$mor) &
  !is.na(network$source) &
  !is.na(network$target),
]

if (anyDuplicated(network[c("source", "target")])) {
  stop("Duplicate TF-target interactions")
}

targets <- function(tf) {
  unique(network$target[network$source == tf])
}

srebf2 <- targets("SREBF2")
srebf1 <- targets("SREBF1")

shared <- intersect(srebf2, srebf1)

sterol_core <- c(
  "TM7SF2", "CYP51A1", "MSMO1", "SQLE",
  "FDFT1", "HSD17B7", "LSS", "NSDHL",
  "DHCR7", "EBP", "SC5D", "DHCR24"
)

measured <- ranked$gene_symbol

cat("=== PHASE 17E: REGULON OVERLAP ===\n")
cat("SREBF2 reference targets:", length(srebf2), "\n")
cat("SREBF1 reference targets:", length(srebf1), "\n")
cat("Shared reference targets:", length(shared), "\n")
cat("Shared measured targets:",
    length(intersect(shared, measured)), "\n")
cat("Shared target names:",
    paste(intersect(shared, measured), collapse = ", "), "\n")

cat("\n=== STEROL-CORE OVERLAP ===\n")
cat("SREBF2 measured core targets:",
    paste(intersect(intersect(srebf2, sterol_core), measured),
          collapse = ", "), "\n")
cat("SREBF1 measured core targets:",
    paste(intersect(intersect(srebf1, sterol_core), measured),
          collapse = ", "), "\n")

# Four scenarios retaining the same DoRothEA reference.
scenarios <- list(
  baseline = character(0),
  exclude_shared_targets = shared,
  exclude_sterol_core = sterol_core,
  exclude_both = union(shared, sterol_core)
)

results <- list()

for (name in names(scenarios)) {
  drop <- scenarios[[name]]

  sub <- ranked[!ranked$gene_symbol %in% drop, ]

  mat <- matrix(
    sub$stat,
    ncol = 1,
    dimnames = list(
      sub$gene_symbol,
      "KO_vs_WT_depletion"
    )
  )

  # Record eligible regulon target counts explicitly.
  net_sub <- network[
    network$target %in% rownames(mat),
  ]

  counts <- table(net_sub$source)

  fit <- decoupleR::run_viper(
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

  fit <- as.data.frame(fit)

  for (tf in c("SREBF2", "SREBF1", "MYC", "SP1")) {
    hit <- fit[fit$source == tf, ]

    n_targets <- if (tf %in% names(counts)) {
      as.integer(counts[[tf]])
    } else {
      0L
    }

    results[[length(results) + 1L]] <- data.frame(
      scenario = name,
      tf = tf,
      remaining_targets = n_targets,
      score = if (nrow(hit)) hit$score[1] else NA_real_,
      p_value = if (nrow(hit)) hit$p_value[1] else NA_real_,
      status = if (nrow(hit)) "ESTIMABLE" else "NOT_ESTIMABLE"
    )
  }
}

out <- do.call(rbind, results)

write.table(
  out,
  file.path(outdir, "phase17e_regulon_dependence.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\n=== PHASE 17E: DEPENDENCE RESULTS ===\n")
print(out, row.names = FALSE, digits = 5)

cat("\nPASS: Phase 17E completed\n")
