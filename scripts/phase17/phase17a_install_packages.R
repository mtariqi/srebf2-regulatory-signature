#!/usr/bin/env Rscript

lib <- path.expand("~/R/srebf2-phase17-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 600
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib)
}

cat("R version:", R.version.string, "\n")
cat("Bioconductor version:",
    as.character(BiocManager::version()), "\n")

cran_packages <- c("ggplot2", "dplyr")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1),
          quietly = TRUE)
]

if (length(missing_cran)) {
  install.packages(missing_cran, lib = lib)
}

bioc_packages <- c(
  "dorothea",
  "viper",
  "decoupleR",
  "BiocParallel"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1),
          quietly = TRUE)
]

if (length(missing_bioc)) {
  BiocManager::install(
    missing_bioc,
    lib = lib,
    ask = FALSE,
    update = FALSE
  )
}

cat("\n=== PHASE 17 DEPENDENCY AUDIT ===\n")

packages <- c(cran_packages, bioc_packages)

for (pkg in packages) {
  installed <- requireNamespace(pkg, quietly = TRUE)

  cat(
    pkg, ":",
    if (installed)
      as.character(packageVersion(pkg))
    else
      "MISSING",
    "\n"
  )
}

if (!all(vapply(
  packages, requireNamespace,
  logical(1), quietly = TRUE
))) {
  stop("Phase 17 dependencies incomplete")
}

cat("\nPASS: Phase 17 dependencies installed\n")
