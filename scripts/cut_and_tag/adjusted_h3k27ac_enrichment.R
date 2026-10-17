#!/usr/bin/env Rscript

# Phase 9: Confounder-adjusted H3K27ac enrichment
# Study: SREBF2 regulatory signature
# Outcome: Presence of H3K27ac peak overlap
# Predictors: Signature, gene length, expression, biotype

input_file <- paste0(
  "results/cut_and_tag/integration/",
  "SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
)

bed_file <- paste0(
  "data/processed/cut_and_tag/annotation/",
  "gencode_v47_gene_bodies.bed"
)

output_dir <- "results/cut_and_tag/integration"
report_dir <- "reports/phase9"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file), file.exists(bed_file))

# ------------------------------------------------------------
# 1. Read RNA-seq / CUT&Tag integration data
# ------------------------------------------------------------

x <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

x <- x[
  x$h3k27ac_match_status %in% c("exact", "stable_id"),
  ,
  drop = FALSE
]

stopifnot(nrow(x) == 17097)
stopifnot(!anyDuplicated(x$gene_id))

# ------------------------------------------------------------
# 2. Read GENCODE gene-body coordinates
# ------------------------------------------------------------

bed <- read.delim(
  bed_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(ncol(bed) >= 7)

colnames(bed)[1:7] <- c(
  "chrom", "start", "end", "gene_id",
  "gene_name", "gene_type", "strand"
)

bed$gene_length <- as.numeric(bed$end) - as.numeric(bed$start)

stopifnot(all(bed$gene_length > 0))
stopifnot(!anyDuplicated(bed$gene_id))

# ------------------------------------------------------------
# 3. Match gene lengths
# ------------------------------------------------------------

idx <- match(x$h3k27ac_gene_id, bed$gene_id)

x$gene_length <- bed$gene_length[idx]

stopifnot(!anyNA(x$gene_length))

# ------------------------------------------------------------
# 4. Prepare model covariates
# ------------------------------------------------------------

stopifnot(all(x$signature_selected %in% c("True", "False")))

x$signature <- as.integer(x$signature_selected == "True")

x$baseMean <- as.numeric(x$baseMean)

stopifnot(all(is.finite(x$baseMean)))
stopifnot(all(x$baseMean >= 0))

x$log_gene_length <- log10(x$gene_length)
x$log_baseMean <- log10(x$baseMean + 1)

x$biotype_group <- ifelse(
  x$gene_type == "protein_coding",
  "protein_coding",
  ifelse(x$gene_type == "lncRNA", "lncRNA", "other")
)

x$biotype_group <- factor(
  x$biotype_group,
  levels = c("protein_coding", "lncRNA", "other")
)

stopifnot(sum(x$signature) == 1316)

# ------------------------------------------------------------
# 5. Fit logistic regression models
# ------------------------------------------------------------

categories <- c(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

results <- list()
models <- list()

for (category in names(categories)) {

  peak_col <- categories[[category]]

  x$overlap <- as.integer(as.numeric(x[[peak_col]]) > 0)

  model <- glm(
    overlap ~ signature +
      log_gene_length +
      log_baseMean +
      biotype_group,
    data = x,
    family = binomial(link = "logit")
  )

  if (!model$converged) {
    stop(paste("Model failed to converge:", category))
  }

  coef_table <- summary(model)$coefficients

  estimate <- coef_table["signature", "Estimate"]
  se <- coef_table["signature", "Std. Error"]
  pvalue <- coef_table["signature", "Pr(>|z|)"]

  results[[category]] <- data.frame(
    category = category,
    n_genes = nrow(x),
    signature_genes = sum(x$signature),
    background_genes = sum(x$signature == 0),
    adjusted_OR = exp(estimate),
    CI_lower = exp(estimate - 1.96 * se),
    CI_upper = exp(estimate + 1.96 * se),
    pvalue = pvalue,
    converged = model$converged
  )

  models[[category]] <- model
}

out <- do.call(rbind, results)
out$padj_BH <- p.adjust(out$pvalue, method = "BH")

out_file <- file.path(
  output_dir,
  "SREBF2_H3K27ac_ADJUSTED_ENRICHMENT.tsv"
)

write.table(
  out,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  models,
  file.path(output_dir, "SREBF2_H3K27ac_ADJUSTED_MODELS.rds")
)

# ------------------------------------------------------------
# 6. QC summary
# ------------------------------------------------------------

qc_file <- file.path(report_dir, "adjusted_model_QC.txt")

sink(qc_file)

cat("PHASE 9: ADJUSTED H3K27ac ENRICHMENT\n\n")
cat("Genes analyzed:", nrow(x), "\n")
cat("Signature genes:", sum(x$signature), "\n")
cat("Non-signature genes:", sum(x$signature == 0), "\n\n")

cat("Biotype distribution:\n")
print(table(x$biotype_group))

cat("\nGene length summary:\n")
print(summary(x$gene_length))

cat("\nBaseline expression summary:\n")
print(summary(x$baseMean))

cat("\nAdjusted enrichment results:\n")
print(out, row.names = FALSE)

sink()

cat("\n========== PHASE 9 RESULTS ==========\n")
print(out, row.names = FALSE)

cat("\nSaved:", out_file, "\n")
cat("QC:", qc_file, "\n")
