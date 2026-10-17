#!/usr/bin/env Rscript

input <- paste0(
  "results/cut_and_tag/integration/",
  "SREBF2_RNASEQ_H3K27ac_BACKGROUND.tsv"
)

output <- paste0(
  "results/cut_and_tag/integration/",
  "SREBF2_H3K27ac_BACKGROUND_ENRICHMENT.tsv"
)

x <- read.delim(
  input,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Restrict to genes with valid annotation matches.
x <- x[
  x$h3k27ac_match_status %in% c("exact", "stable_id"),
]

stopifnot(nrow(x) == 17097)

selected <- x$signature_selected == "True"

stopifnot(sum(selected) == 1316)
stopifnot(sum(!selected) == 15781)

categories <- c(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

results <- list()

for (category in names(categories)) {

  overlap <- as.integer(x[[categories[[category]]]]) > 0

  tab <- matrix(
    c(
      sum(selected & overlap),
      sum(selected & !overlap),
      sum(!selected & overlap),
      sum(!selected & !overlap)
    ),
    nrow = 2,
    byrow = TRUE
  )

  test <- fisher.test(tab, alternative = "two.sided")

  results[[category]] <- data.frame(
    category = category,
    selected_overlap = tab[1, 1],
    selected_no_overlap = tab[1, 2],
    nonselected_overlap = tab[2, 1],
    nonselected_no_overlap = tab[2, 2],
    odds_ratio = unname(test$estimate),
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    p_value = test$p.value
  )
}

out <- do.call(rbind, results)
out$BH_FDR <- p.adjust(out$p_value, method = "BH")

write.table(
  out,
  output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(out, row.names = FALSE)

cat("\nSaved:", output, "\n")
