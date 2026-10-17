#!/usr/bin/env Rscript

# Phase 10: SREBF2 / H3K27ac enrichment sensitivity analyses

input <- "results/cut_and_tag/integration/SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
bed_file <- "data/processed/cut_and_tag/annotation/gencode_v47_gene_bodies.bed"
output <- "results/cut_and_tag/integration/SREBF2_H3K27ac_PHASE10_SENSITIVITY.tsv"

stopifnot(file.exists(input), file.exists(bed_file))

x <- read.delim(input, stringsAsFactors = FALSE, check.names = FALSE)

x <- x[x$h3k27ac_match_status %in% c("exact", "stable_id"), , drop = FALSE]

stopifnot(nrow(x) == 17097, !anyDuplicated(x$gene_id))

bed <- read.delim(bed_file, header = FALSE, stringsAsFactors = FALSE)

stopifnot(ncol(bed) >= 7)

names(bed)[1:7] <- c(
  "chrom", "start", "end", "gene_id",
  "gene_name", "gene_type", "strand"
)

stopifnot(!anyDuplicated(bed$gene_id))

bed$gene_length <- as.numeric(bed$end) - as.numeric(bed$start)

x$gene_length <- bed$gene_length[
  match(x$h3k27ac_gene_id, bed$gene_id)
]

x$baseMean <- as.numeric(x$baseMean)

stopifnot(
  all(is.finite(x$gene_length)),
  all(x$gene_length > 0),
  all(is.finite(x$baseMean)),
  all(x$baseMean >= 0),
  all(x$signature_selected %in% c("True", "False"))
)

x$signature <- as.integer(x$signature_selected == "True")
x$log_gene_length <- log10(x$gene_length)
x$log_baseMean <- log10(x$baseMean + 1)

x$biotype_group <- factor(
  ifelse(
    x$gene_type == "protein_coding",
    "protein_coding",
    ifelse(x$gene_type == "lncRNA", "lncRNA", "other")
  ),
  levels = c("protein_coding", "lncRNA", "other")
)

categories <- c(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

results <- list()

fit_model <- function(dat, category, analysis, formula) {

  dat$overlap <- as.integer(
    as.numeric(dat[[categories[[category]]]]) > 0
  )

  if (length(unique(dat$signature)) != 2) {
    stop(paste("Both signature groups required:", analysis))
  }

  model <- glm(formula, data = dat, family = binomial())

  if (!model$converged) {
    stop(paste("Nonconverged model:", analysis, category))
  }

  cf <- summary(model)$coefficients

  if (!("signature" %in% rownames(cf))) {
    stop("Signature coefficient missing")
  }

  b <- cf["signature", "Estimate"]
  se <- cf["signature", "Std. Error"]

  data.frame(
    analysis = analysis,
    category = category,
    n_genes = nrow(dat),
    n_signature = sum(dat$signature),
    n_background = sum(dat$signature == 0),
    adjusted_OR = exp(b),
    CI_lower = exp(b - 1.96 * se),
    CI_upper = exp(b + 1.96 * se),
    pvalue = cf["signature", "Pr(>|z|)"],
    converged = model$converged
  )
}

# Sensitivity A: protein-coding genes only

pc <- x[x$gene_type == "protein_coding", , drop = FALSE]

for (category in names(categories)) {
  results[[length(results) + 1]] <- fit_model(
    pc,
    category,
    "protein_coding_only",
    overlap ~ signature + log_gene_length + log_baseMean
  )
}

# Sensitivity B: expression quartiles
# Quartile boundaries are calculated across all matched genes.

x$expression_quartile <- cut(
  rank(x$baseMean, ties.method = "first"),
  breaks = 4,
  labels = paste0("Q", 1:4)
)

for (q in levels(x$expression_quartile)) {

  subset <- x[x$expression_quartile == q, , drop = FALSE]
  subset$biotype_group <- droplevels(subset$biotype_group)

  for (category in names(categories)) {
    results[[length(results) + 1]] <- fit_model(
      subset,
      category,
      paste0("expression_", q),
      overlap ~ signature + log_gene_length + biotype_group
    )
  }
}

# Sensitivity C: flexible gene-length adjustment
# Natural cubic splines allow a nonlinear length-overlap relationship.

for (category in names(categories)) {
  results[[length(results) + 1]] <- fit_model(
    x,
    category,
    "gene_length_spline",
    overlap ~ signature +
      splines::ns(log_gene_length, df = 3) +
      log_baseMean + biotype_group
  )
}

out <- do.call(rbind, results)

out$padj_BH <- p.adjust(out$pvalue, method = "BH")

write.table(
  out,
  output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========== PHASE 10 SENSITIVITY RESULTS ==========\n")
print(out, row.names = FALSE)

cat("\nSaved:", output, "\n")
