#!/usr/bin/env Rscript

lib <- path.expand("~/R/srebf2-phase17-library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(dorothea)
  library(decoupleR)
  library(ggplot2)
})

outdir <- "results/phase17"

ranked <- read.delim(
  file.path(outdir, "phase17b_ranked_unique_symbols.tsv")
)

data("dorothea_hs", package = "dorothea")

reg <- as.data.frame(dorothea_hs)
reg <- reg[reg$confidence %in% c("A", "B", "C"), ]

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

stopifnot(!anyDuplicated(network[c("source", "target")]))

shared <- intersect(
  network$target[network$source == "SREBF2"],
  network$target[network$source == "SREBF1"]
)

core <- c(
  "TM7SF2", "CYP51A1", "MSMO1", "SQLE",
  "FDFT1", "HSD17B7", "LSS", "NSDHL",
  "DHCR7", "EBP", "SC5D", "DHCR24"
)

scenarios <- list(
  Baseline = character(0),
  "Exclude shared" = shared,
  "Exclude sterol core" = core,
  "Exclude both" = union(shared, core)
)

records <- list()

for (name in names(scenarios)) {
  sub <- ranked[
    !ranked$gene_symbol %in% scenarios[[name]],
  ]

  mat <- matrix(
    sub$stat,
    ncol = 1,
    dimnames = list(
      sub$gene_symbol,
      "KO_vs_WT_depletion"
    )
  )

  net <- network[
    network$target %in% rownames(mat),
  ]

  counts <- table(net$source)

  # Five-target exploratory sensitivity analysis
  fit <- as.data.frame(decoupleR::run_viper(
    mat = mat,
    network = network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = 5,
    method = "none",
    pleiotropy = FALSE,
    verbose = FALSE
  ))

  for (tf in c("SREBF2", "SREBF1")) {
    hit <- fit[fit$source == tf, ]

    n <- if (tf %in% names(counts)) {
      as.integer(counts[[tf]])
    } else 0L

    records[[length(records) + 1L]] <- data.frame(
      scenario = name,
      tf = tf,
      measured_targets = n,
      score = if (nrow(hit)) hit$score[1] else NA_real_,
      p_value = if (nrow(hit)) hit$p_value[1] else NA_real_,
      status = if (nrow(hit)) "ESTIMABLE" else "NOT_ESTIMABLE"
    )
  }
}

out <- do.call(rbind, records)

out$scenario <- factor(
  out$scenario,
  levels = names(scenarios)
)

write.table(
  out,
  file.path(outdir, "phase17f_five_target_sensitivity.tsv"),
  sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)

cat("=== PHASE 17F: FIVE-TARGET SENSITIVITY ===\n")
print(out, row.names = FALSE)

p <- ggplot(
  out[!is.na(out$score), ],
  aes(x = scenario, y = score, fill = tf)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  facet_wrap(~tf, ncol = 1) +
  labs(
    title = "SREBF regulatory activity sensitivity",
    subtitle = paste(
      "DoRothEA A-C; VIPER; minimum five targets",
      "| KO vs WT under cholesterol depletion"
    ),
    x = NULL,
    y = "Contrast-level VIPER score",
    caption = paste(
      "Missing bars indicate insufficient target coverage;",
      "scores are model-based, not experimental replication."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 20, hjust = 1
    ),
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "Figure_TF_Activity_Sensitivity.pdf"),
  p, width = 8.5, height = 6
)

ggsave(
  file.path(outdir, "Figure_TF_Activity_Sensitivity.png"),
  p, width = 8.5, height = 6, dpi = 300
)

cat("\nPASS: Phase 17F completed\n")
