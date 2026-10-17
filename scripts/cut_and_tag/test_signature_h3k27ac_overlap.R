#!/usr/bin/env Rscript

input <- "results/cut_and_tag/integration/SREBF2_SIGNATURE_H3K27ac_INTEGRATED.tsv"

output <- "results/cut_and_tag/integration/SREBF2_H3K27ac_Fisher_results.tsv"

x <- read.delim(
  input,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

x <- x[
  x$h3k27ac_match_status %in% c("exact", "stable_id"),
]

stopifnot(nrow(x) == 1316)

groups <- c("increased_in_KO", "decreased_in_KO")

x <- x[x$direction %in% groups, ]

tests <- list(
  promoter = "promoter_peak_count",
  gene_body = "gene_body_peak_count",
  any = "any_peak_count"
)

results <- list()

for (name in names(tests)) {

  column <- tests[[name]]

  overlap <- as.integer(x[[column]]) > 0

  tab <- table(
    factor(x$direction, levels = groups),
    factor(overlap, levels = c(TRUE, FALSE))
  )

  result <- fisher.test(tab)

  results[[name]] <- data.frame(
    comparison = name,
    increased_overlap = tab[1, 1],
    increased_no_overlap = tab[1, 2],
    decreased_overlap = tab[2, 1],
    decreased_no_overlap = tab[2, 2],
    odds_ratio = unname(result$estimate),
    ci_lower = result$conf.int[1],
    ci_upper = result$conf.int[2],
    p_value = result$p.value
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
