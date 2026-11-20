#!/usr/bin/env Rscript

lib <- path.expand("~/R/srebf2-phase17-library")
stopifnot(dir.exists(lib))
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

mat <- matrix(
  ranked$stat,
  ncol = 1,
  dimnames = list(
    ranked$gene_symbol,
    "KO_vs_WT_depletion"
  )
)

data("dorothea_hs", package = "dorothea")

reg <- as.data.frame(dorothea_hs)

configs <- list(
  AC_10 = list(conf = c("A", "B", "C"), minsize = 10),
  AC_5  = list(conf = c("A", "B", "C"), minsize = 5),
  AB_10 = list(conf = c("A", "B"), minsize = 10),
  AB_5  = list(conf = c("A", "B"), minsize = 5)
)

all_results <- list()

for (name in names(configs)) {
  cfg <- configs[[name]]

  sub <- reg[
    reg$confidence %in% cfg$conf &
    reg$target %in% rownames(mat),
  ]

  net <- data.frame(
    source = sub$tf,
    target = sub$target,
    mor = as.numeric(sub$mor)
  )

  net <- net[
    is.finite(net$mor) &
    !is.na(net$source) &
    !is.na(net$target),
  ]

  if (anyDuplicated(net[c("source", "target")])) {
    stop(paste("Duplicate TF-target edges:", name))
  }

  if (nrow(net) == 0) {
    stop(paste("Empty regulon:", name))
  }

  set.seed(1704)

  fit <- decoupleR::run_viper(
    mat = mat,
    network = net,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = cfg$minsize,
    method = "none",
    pleiotropy = FALSE,
    verbose = FALSE
  )

  fit <- as.data.frame(fit)

  if (!nrow(fit)) {
    stop(paste("No TF activity estimates:", name))
  }

  if (anyDuplicated(fit$source)) {
    stop(paste("Duplicate TF results:", name))
  }

  fit$fdr_bh <- p.adjust(
    fit$p_value,
    method = "BH"
  )

  target_counts <- aggregate(
    target ~ source,
    data = unique(net[c("source", "target")]),
    FUN = length
  )

  names(target_counts)[2] <- "measured_targets"

  fit <- merge(
    fit,
    target_counts,
    by = "source",
    all.x = TRUE
  )

  fit$configuration <- name
  fit$confidence <- paste(cfg$conf, collapse = "")
  fit$minsize <- cfg$minsize

  all_results[[name]] <- fit

  cat(
    name,
    ": TFs=", nrow(fit),
    " interactions=", nrow(net),
    "\n",
    sep = ""
  )
}

combined <- do.call(rbind, all_results)

write.table(
  combined,
  file.path(outdir, "phase17d_viper_sensitivity.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

candidates <- c(
  "SREBF2", "SREBF1", "MYC", "SP1",
  "ZBTB7A", "ZNF263", "MAZ", "IRF3"
)

grid <- expand.grid(
  source = candidates,
  configuration = names(configs),
  stringsAsFactors = FALSE
)

candidate_results <- merge(
  grid,
  combined[, c(
    "source", "configuration",
    "score", "p_value",
    "fdr_bh", "measured_targets"
  )],
  by = c("source", "configuration"),
  all.x = TRUE,
  sort = FALSE
)

candidate_results <- candidate_results[
  order(
    match(candidate_results$source, candidates),
    match(candidate_results$configuration, names(configs))
  ),
]

candidate_results$direction <- ifelse(
  is.na(candidate_results$score),
  "NOT_ESTIMABLE",
  ifelse(
    candidate_results$score < 0,
    "NEGATIVE",
    ifelse(candidate_results$score > 0, "POSITIVE", "ZERO")
  )
)

write.table(
  candidate_results,
  file.path(outdir, "phase17d_candidate_tf_sensitivity.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

cat("\n=== PHASE 17D: CANDIDATE TF SENSITIVITY ===\n")

print(
  candidate_results,
  row.names = FALSE,
  digits = 5
)

cat("\n=== PRIMARY A-C, MINIMUM 10: REPRODUCIBILITY ===\n")

original <- read.delim(
  file.path(outdir, "phase17c_viper_activity.tsv")
)

primary <- all_results[["AC_10"]]

check <- merge(
  original[, c("source", "score")],
  primary[, c("source", "score")],
  by = "source",
  suffixes = c("_original", "_rerun")
)

stopifnot(
  nrow(check) == nrow(original),
  nrow(check) == nrow(primary)
)

max_diff <- max(abs(
  check$score_original - check$score_rerun
))

cat("Maximum absolute score difference:", max_diff, "\n")

if (!isTRUE(all.equal(
  check$score_original,
  check$score_rerun,
  tolerance = 1e-8
))) {
  stop("Primary VIPER scores did not reproduce")
}

cat("\nPASS: Phase 17D complete\n")
