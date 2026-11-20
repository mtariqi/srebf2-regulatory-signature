#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)


# Phase 17 isolated R package library
phase17_lib <- path.expand("~/R/srebf2-phase17-library")

if (!dir.exists(phase17_lib)) {
    stop("Phase 17 R library directory missing")
}

.libPaths(c(phase17_lib, .libPaths()))

library(dorothea)

dir.create(
  "results/phase17",
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. Load verified DESeq2 Wald-ranked data
ranked <- read.delim(
  "results/phase15/gsea_primary_wald_ranked.tsv",
  check.names = FALSE
)

stopifnot(
  nrow(ranked) == 12485,
  all(c("gene_id", "gene_symbol", "stat") %in% names(ranked)),
  !anyDuplicated(ranked$gene_id)
)

ranked$stat <- as.numeric(ranked$stat)

ranked <- ranked[
  is.finite(ranked$stat) &
  !is.na(ranked$gene_symbol) &
  nzchar(trimws(ranked$gene_symbol)),
]

# 2. Exclude ambiguously mapped gene symbols
symbol_counts <- table(ranked$gene_symbol)

duplicate_symbols <- names(
  symbol_counts[symbol_counts > 1]
)

cat("=== GENE IDENTIFIER AUDIT ===\n")
cat("Valid ranked genes:", nrow(ranked), "\n")
cat("Ambiguous gene symbols:", length(duplicate_symbols), "\n")

ranked <- ranked[
  !ranked$gene_symbol %in% duplicate_symbols,
]

stopifnot(!anyDuplicated(ranked$gene_symbol))

# 3. Load human DoRothEA regulons
data("dorothea_hs", package = "dorothea")

regulons <- dorothea_hs

stopifnot(
  all(c("tf", "target", "confidence") %in% names(regulons))
)

regulons <- regulons[
  regulons$confidence %in% c("A", "B", "C"),
]

stopifnot(nrow(regulons) > 0)

# 4. Compute target coverage per TF
tf_total <- aggregate(
  target ~ tf,
  data = unique(regulons[c("tf", "target")]),
  FUN = length
)

names(tf_total)[2] <- "total_targets"

overlap <- regulons[
  regulons$target %in% ranked$gene_symbol,
]

tf_measured <- aggregate(
  target ~ tf,
  data = unique(overlap[c("tf", "target")]),
  FUN = length
)

names(tf_measured)[2] <- "measured_targets"

coverage <- merge(
  tf_total,
  tf_measured,
  by = "tf",
  all.x = TRUE
)

coverage$measured_targets[
  is.na(coverage$measured_targets)
] <- 0L

coverage$coverage_fraction <- (
  coverage$measured_targets /
  coverage$total_targets
)

coverage$eligible_5_targets <- (
  coverage$measured_targets >= 5
)

coverage$eligible_10_targets <- (
  coverage$measured_targets >= 10
)

coverage <- coverage[
  order(-coverage$measured_targets, coverage$tf),
]

# 5. Save auditable outputs
write.table(
  coverage,
  "results/phase17/phase17b_regulon_coverage.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  ranked,
  "results/phase17/phase17b_ranked_unique_symbols.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 6. Summaries
cat("\n=== DOROTHEA HUMAN REGULON AUDIT ===\n")
cat("A-C regulatory interactions:", nrow(regulons), "\n")
cat("Unique transcription factors:", length(unique(regulons$tf)), "\n")
cat("TFs with >= 5 measured targets:",
    sum(coverage$eligible_5_targets), "\n")
cat("TFs with >= 10 measured targets:",
    sum(coverage$eligible_10_targets), "\n")

cat("\n=== CANDIDATE REGULATORS ===\n")

candidates <- c("SREBF2", "SREBF1", "MYC", "SP1")

for (tf in candidates) {
  x <- coverage[coverage$tf == tf, ]

  if (nrow(x) == 0) {
    cat(tf, ": NOT REPRESENTED\n")
  } else {
    cat(
      tf,
      ": measured=", x$measured_targets,
      " total=", x$total_targets,
      " coverage=", round(x$coverage_fraction, 3),
      "\n",
      sep = ""
    )
  }
}

cat("\n=== TOP 15 TFs BY TARGET COVERAGE ===\n")
print(head(coverage, 15), row.names = FALSE)

cat("\nPASS: Phase 17B coverage audit complete\n")
