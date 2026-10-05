conda <- Sys.getenv("CONDA_PREFIX")
stopifnot(nzchar(conda))
.libPaths(file.path(conda, "lib/R/library"), include.site = FALSE)
suppressPackageStartupMessages(library(DESeq2))

out <- "results/qc/validation/GSE271000/expression"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(
  "results/counts/validation/GSE271000/GSE271000_validation_gene_counts.tsv",
  row.names = 1, check.names = FALSE
)
samples <- read.delim(
  "metadata/GSE271000_validation_samples.tsv",
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(counts) == 78932L,
  ncol(counts) == 8L,
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(samples$run),
  identical(colnames(counts), samples$run)
)

counts <- as.matrix(counts)
stopifnot(
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == floor(counts))
)
storage.mode(counts) <- "integer"
rownames(samples) <- samples$run

# Expression filter for cohort QC only.
keep <- rowSums(counts >= 10L) >= 2L
stopifnot(sum(keep) > 500L)

dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = samples,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
x <- assay(vsd)

write_tsv <- function(data, name, row_names = FALSE) {
  write.table(
    data, file.path(out, name), sep = "\t",
    quote = FALSE, row.names = row_names,
    col.names = if (row_names) NA else TRUE
  )
}

summary <- data.frame(
  run = samples$run,
  group = samples$group,
  assigned_counts = colSums(counts),
  nonzero_genes = colSums(counts > 0),
  size_factor = sizeFactors(dds)
)
write_tsv(summary, "SAMPLE_NORMALIZATION_QC.tsv")
write_tsv(
  data.frame(gene_id = rownames(counts), retained_for_QC = keep),
  "QC_GENE_FILTER.tsv"
)

# PCA uses the 500 most variable retained genes.
variances <- apply(x, 1, var)
selected <- head(order(variances, decreasing = TRUE), 500L)
pca <- prcomp(t(x[selected, , drop = FALSE]), center = TRUE)
percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(
  run = samples$run,
  group = samples$group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
write_tsv(scores, "PCA_SCORES.tsv")
write_tsv(
  data.frame(PC = seq_along(percent), variance_percent = percent),
  "PCA_VARIANCE.tsv"
)

# Correlations use all retained genes after transformation.
correlations <- cor(x, method = "pearson")
write_tsv(as.data.frame(correlations), "SAMPLE_CORRELATIONS.tsv", TRUE)

pairs <- combn(seq_len(ncol(x)), 2L)
pair_table <- data.frame(
  run_1 = samples$run[pairs[1, ]],
  run_2 = samples$run[pairs[2, ]],
  same_group = samples$group[pairs[1, ]] == samples$group[pairs[2, ]],
  correlation = correlations[t(pairs)]
)
write_tsv(pair_table, "PAIRWISE_CORRELATIONS.tsv")

groups <- c("control", "shSREBP2_664", "shSREBP2_665", "shSREBP2_667")
palette <- setNames(c("#333333", "#0072B2", "#D55E00", "#009E73"), groups)
colors <- unname(palette[samples$group])
labels <- paste(samples$group, samples$replicate_label, sep = "_")

pdf(file.path(out, "VALIDATION_EXPRESSION_QC.pdf"), width = 10, height = 8)

barplot(
  colSums(counts) / 1e6, names.arg = labels, col = colors,
  las = 2, ylab = "Assigned counts (millions)",
  main = "Validation library sizes", cex.names = 0.7
)

plot(
  scores$PC1, scores$PC2, col = colors, pch = 19, cex = 1.5,
  xlab = sprintf("PC1 (%.2f%%)", percent[1]),
  ylab = sprintf("PC2 (%.2f%%)", percent[2]),
  main = "Blind VST PCA: top 500 variable genes"
)
text(scores$PC1, scores$PC2, labels = labels, pos = 3, cex = 0.65)
legend("topright", legend = groups, col = palette, pch = 19, cex = 0.8)

heatmap(
  correlations, symm = TRUE, scale = "none",
  labRow = labels, labCol = labels,
  margins = c(10, 10), cexRow = 0.8, cexCol = 0.8,
  main = "Sample correlations: all QC-retained genes"
)
dev.off()

writeLines(
  trimws(capture.output(sessionInfo()), which = "right"),
  file.path(out, "R_SESSION_INFO.txt")
)
writeLines(c(
  "Purpose: cohort expression QC; no differential-expression tests.",
  "QC filter: count >=10 in at least two of eight samples.",
  "Normalization: DESeq2 median-ratio size factors.",
  "Transformation: varianceStabilizingTransformation, blind=TRUE.",
  "PCA: centered transformed expression of top 500 variable genes.",
  "Correlations: Pearson, all QC-retained transformed genes.",
  "No samples automatically excluded.",
  "Signature definition must be frozen before validation-effect analysis."
), file.path(out, "METHODS.txt"))

cat("Input genes:", nrow(counts), "\n")
cat("Genes retained for QC:", sum(keep), "\n")
cat("PC1 variance:", round(percent[1], 2), "%\n")
cat("PC2 variance:", round(percent[2], 2), "%\n")
print(summary, row.names = FALSE)
cat("\nWithin-group correlations:\n")
print(pair_table[pair_table$same_group, ], row.names = FALSE)
cat("\nSaved:", out, "\n")
