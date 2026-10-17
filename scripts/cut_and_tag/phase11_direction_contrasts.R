#!/usr/bin/env Rscript

input <- "results/cut_and_tag/integration/SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
bed_file <- "data/processed/cut_and_tag/annotation/gencode_v47_gene_bodies.bed"
output <- "results/cut_and_tag/integration/SREBF2_H3K27ac_DIRECTION_CONTRASTS.tsv"

x <- read.delim(input, stringsAsFactors = FALSE)
bed <- read.delim(bed_file, header = FALSE)

names(bed)[1:7] <- c(
  "chrom", "start", "end", "gene_id",
  "gene_name", "gene_type", "strand"
)

x <- x[x$h3k27ac_match_status %in% c("exact", "stable_id"), ]

stopifnot(
  nrow(x) == 17097,
  !anyDuplicated(x$gene_id),
  !anyDuplicated(bed$gene_id)
)

lengths <- as.numeric(bed$end) - as.numeric(bed$start)

x$gene_length <- lengths[
  match(x$h3k27ac_gene_id, bed$gene_id)
]

x$baseMean <- as.numeric(x$baseMean)

stopifnot(
  all(is.finite(x$gene_length)),
  all(x$gene_length > 0),
  all(is.finite(x$baseMean)),
  all(x$baseMean >= 0)
)

x$log_gene_length <- log10(x$gene_length)
x$log_baseMean <- log10(x$baseMean + 1)

x$biotype_group <- factor(
  ifelse(
    x$gene_type == "protein_coding",
    "protein_coding",
    ifelse(x$gene_type == "lncRNA", "lncRNA", "other")
  )
)

x$direction_group <- factor(
  x$direction,
  levels = c(
    "not_selected",
    "decreased_in_KO",
    "increased_in_KO"
  )
)

categories <- c(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

results <- list()

for (category in names(categories)) {

  x$overlap <- as.integer(
    as.numeric(x[[categories[[category]]]]) > 0
  )

  model <- glm(
    overlap ~ direction_group +
      log_gene_length +
      log_baseMean +
      biotype_group,
    data = x,
    family = binomial()
  )

  stopifnot(model$converged)

  beta <- coef(model)
  covariance <- vcov(model)

  a <- "direction_groupincreased_in_KO"
  b <- "direction_groupdecreased_in_KO"

  difference <- beta[a] - beta[b]

  variance <- covariance[a, a] +
    covariance[b, b] -
    2 * covariance[a, b]

  stopifnot(is.finite(variance), variance > 0)

  se <- sqrt(variance)
  z <- difference / se
  pvalue <- 2 * pnorm(-abs(z))

  results[[category]] <- data.frame(
    category = category,
    comparison = "increased_vs_decreased",
    adjusted_OR = unname(exp(difference)),
    CI_lower = unname(exp(difference - 1.96 * se)),
    CI_upper = unname(exp(difference + 1.96 * se)),
    pvalue = unname(pvalue),
    converged = model$converged
  )
}

out <- do.call(rbind, results)
out$padj_BH <- p.adjust(out$pvalue, method = "BH")

write.table(
  out,
  output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n========== DIRECT DIRECTION CONTRASTS ==========\n")
print(out, row.names = FALSE)

cat("\nSaved:", output, "\n")
