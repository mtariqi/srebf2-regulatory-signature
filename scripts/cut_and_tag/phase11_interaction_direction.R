#!/usr/bin/env Rscript

# Phase 11: Expression interaction and direction-specific
# H3K27ac enrichment in the SREBF2 signature.

input <- paste0(
  "results/cut_and_tag/integration/",
  "SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
)

bed_file <- paste0(
  "data/processed/cut_and_tag/annotation/",
  "gencode_v47_gene_bodies.bed"
)

outdir <- "results/cut_and_tag/integration"
reportdir <- "reports/phase11"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reportdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(input, stringsAsFactors = FALSE)
bed <- read.delim(bed_file, header = FALSE)

names(bed)[1:7] <- c(
  "chrom", "start", "end", "gene_id",
  "gene_name", "gene_type", "strand"
)

stopifnot(!anyDuplicated(bed$gene_id))

x <- x[
  x$h3k27ac_match_status %in% c("exact", "stable_id"),
  ,
  drop = FALSE
]

stopifnot(nrow(x) == 17097, !anyDuplicated(x$gene_id))

x$gene_length <- (
  bed$end - bed$start
)[match(x$h3k27ac_gene_id, bed$gene_id)]

x$baseMean <- as.numeric(x$baseMean)

stopifnot(
  all(is.finite(x$gene_length)),
  all(x$gene_length > 0),
  all(is.finite(x$baseMean)),
  all(x$baseMean >= 0)
)

x$log_gene_length <- log10(x$gene_length)
x$log_baseMean <- log10(x$baseMean + 1)

x$signature <- as.integer(x$signature_selected == "True")

x$biotype_group <- factor(
  ifelse(
    x$gene_type == "protein_coding",
    "protein_coding",
    ifelse(x$gene_type == "lncRNA", "lncRNA", "other")
  )
)

x$expression_quartile <- cut(
  rank(x$baseMean, ties.method = "first"),
  breaks = 4,
  labels = paste0("Q", 1:4)
)

x$direction_group <- factor(
  x$direction,
  levels = c(
    "not_selected",
    "decreased_in_KO",
    "increased_in_KO"
  )
)

stopifnot(
  sum(x$direction_group == "not_selected") == 15781,
  sum(x$direction_group == "decreased_in_KO") == 789,
  sum(x$direction_group == "increased_in_KO") == 527
)

categories <- c(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

interaction_results <- list()
direction_results <- list()

for (category in names(categories)) {

  x$overlap <- as.integer(
    as.numeric(x[[categories[[category]]]]) > 0
  )

  # A. Formal signature-by-expression interaction

  reduced <- glm(
    overlap ~ signature + expression_quartile +
      log_gene_length + biotype_group,
    data = x,
    family = binomial()
  )

  full <- glm(
    overlap ~ signature * expression_quartile +
      log_gene_length + biotype_group,
    data = x,
    family = binomial()
  )

  stopifnot(reduced$converged, full$converged)

  comparison <- anova(
    reduced,
    full,
    test = "Chisq"
  )

  interaction_results[[category]] <- data.frame(
    category = category,
    n_genes = nrow(x),
    LR_chisq = comparison$Deviance[2],
    df = comparison$Df[2],
    pvalue = comparison$`Pr(>Chi)`[2],
    reduced_AIC = AIC(reduced),
    full_AIC = AIC(full),
    converged = TRUE
  )

  # B. Direction-specific adjusted enrichment

  direction_model <- glm(
    overlap ~ direction_group +
      log_gene_length +
      log_baseMean +
      biotype_group,
    data = x,
    family = binomial()
  )

  stopifnot(direction_model$converged)

  cf <- summary(direction_model)$coefficients

  for (direction in c(
    "decreased_in_KO",
    "increased_in_KO"
  )) {

    term <- paste0("direction_group", direction)

    b <- cf[term, "Estimate"]
    se <- cf[term, "Std. Error"]

    direction_results[[length(direction_results) + 1]] <-
      data.frame(
        category = category,
        direction = direction,
        n_genes = nrow(x),
        n_direction = sum(x$direction_group == direction),
        adjusted_OR = exp(b),
        CI_lower = exp(b - 1.96 * se),
        CI_upper = exp(b + 1.96 * se),
        pvalue = cf[term, "Pr(>|z|)"],
        converged = direction_model$converged
      )
  }
}

interaction_out <- do.call(rbind, interaction_results)
direction_out <- do.call(rbind, direction_results)

interaction_out$padj_BH <- p.adjust(
  interaction_out$pvalue,
  method = "BH"
)

# Six direction-specific coefficient tests across
# three outcomes and two signature directions.
direction_out$padj_BH <- p.adjust(
  direction_out$pvalue,
  method = "BH"
)

interaction_file <- file.path(
  outdir,
  "SREBF2_H3K27ac_EXPRESSION_INTERACTION.tsv"
)

direction_file <- file.path(
  outdir,
  "SREBF2_H3K27ac_DIRECTION_ADJUSTED.tsv"
)

write.table(
  interaction_out,
  interaction_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  direction_out,
  direction_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n========== EXPRESSION INTERACTION ==========\n")
print(interaction_out, row.names = FALSE)

cat("\n========== DIRECTION-SPECIFIC ENRICHMENT ==========\n")
print(direction_out, row.names = FALSE)

cat("\nSaved:\n", interaction_file, "\n", direction_file, "\n")
