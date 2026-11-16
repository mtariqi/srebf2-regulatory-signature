#!/usr/bin/env Rscript
# ==============================================================================
# Project     : Chromatin-state context of SREBF2 knockout-responsive genes in HeLa cells: an exploratory integrative epigenomic analysis
# File        : scripts/phase15/phase15e_enrichment.R
# Purpose     : GO/Reactome/KEGG ORA and preranked GSEA from frozen Phase 15C
# Author      : Md Tariqul Islam
# Institution : Northeastern University | Bioinformatics Department
# License     : MIT (see LICENSE)
# ==============================================================================
options(stringsAsFactors = FALSE)
set.seed(1515)
required <- c("clusterProfiler", "ReactomePA", "org.Hs.eg.db", "AnnotationDbi", "BiocManager")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
root <- "results/phase15"
out <- file.path(root, "phase15e_enrichment")
if (file.exists(out)) stop("Refusing to overwrite existing results: ", out)
files <- file.path(root, c("protein_coding_universe.tsv", "increased_in_KO.tsv",
                           "decreased_in_KO.tsv", "gsea_primary_wald_ranked.tsv"))
names(files) <- c("universe", "increased", "decreased", "ranked")
if (!all(file.exists(files))) stop("Missing input(s): ", paste(files[!file.exists(files)], collapse = ", "))
dat <- lapply(files, function(p) read.delim(p, check.names = FALSE, quote = "", comment.char = ""))
getcol <- function(d, candidates) {
  x <- intersect(candidates, names(d))
  if (!length(x)) stop("Missing expected column. Found: ", paste(names(d), collapse = ", "))
  x[1]
}
ids <- lapply(dat, function(d) {
  x <- as.character(d[[getcol(d, c("gene_id", "ensembl_gene_id", "ENSEMBL"))]])
  if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) stop("Invalid or duplicate gene identifiers")
  x
})
stopifnot(length(ids$universe) == 12485L, length(ids$increased) == 330L,
          length(ids$decreased) == 625L, length(ids$ranked) == 12485L,
          setequal(ids$universe, ids$ranked),
          all(ids$increased %in% ids$universe),
          all(ids$decreased %in% ids$universe),
          !any(ids$increased %in% ids$decreased))
statcol <- getcol(dat$ranked, c("stat", "wald_stat", "Wald_stat", "wald_statistic", "DESeq2_Wald_stat", "rank_metric"))
wald <- suppressWarnings(as.numeric(dat$ranked[[statcol]]))
if (any(!is.finite(wald))) stop("Non-finite DESeq2 Wald statistics")
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
ens <- strip_version(ids$universe)
if (anyDuplicated(ens)) stop("Duplicate Ensembl IDs after version stripping")
raw <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ens,
                              columns = "ENTREZID", keytype = "ENSEMBL")
raw <- unique(raw[!is.na(raw$ENTREZID) & nzchar(raw$ENTREZID), c("ENSEMBL", "ENTREZID")])
ambiguous_ens <- names(which(table(raw$ENSEMBL) > 1L))
map <- raw[!raw$ENSEMBL %in% ambiguous_ens, , drop = FALSE]
colliding_entrez <- names(which(table(map$ENTREZID) > 1L))
map <- map[!map$ENTREZID %in% colliding_entrez, , drop = FALSE]
if (nrow(map) < 5000 || anyDuplicated(map$ENSEMBL) || anyDuplicated(map$ENTREZID))
  stop("Identifier mapping failed quality checks")
mapids <- function(x) {
  v <- map$ENTREZID[match(strip_version(x), map$ENSEMBL)]
  unique(v[!is.na(v)])
}
universe <- mapids(ids$universe)
up <- mapids(ids$increased)
down <- mapids(ids$decreased)
index <- match(strip_version(ids$ranked), map$ENSEMBL)
valid <- !is.na(index)
rankvec <- wald[valid]
names(rankvec) <- map$ENTREZID[index[valid]]
rankvec <- sort(rankvec, decreasing = TRUE)
stopifnot(!anyDuplicated(names(rankvec)), setequal(names(rankvec), universe))
dir.create(out, recursive = TRUE)
write.table(map, file.path(out, "ensembl_entrez_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(ENTREZID = names(rankvec), Wald_stat = unname(rankvec)),
            file.path(out, "ranked_entrez_wald.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
qa <- data.frame(metric = c("universe_total", "universe_mapped", "increased_total",
                            "increased_mapped", "decreased_total", "decreased_mapped",
                            "ranked_total", "ranked_mapped", "ambiguous_ensembl",
                            "colliding_entrez"),
                 count = c(length(ens), length(universe), length(ids$increased),
                           length(up), length(ids$decreased), length(down),
                           length(ids$ranked), length(rankvec),
                           length(ambiguous_ens), length(colliding_entrez)))
write.table(qa, file.path(out, "mapping_qa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(qa, row.names = FALSE)
run_one <- function(label, fn) {
  cat("RUN:", label, "\n")
  tryCatch({
    result <- fn()
    if (is.null(result)) stop("No result returned")
    tab <- as.data.frame(result)
    write.table(tab, file.path(out, paste0(label, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    saveRDS(result, file.path(out, paste0(label, ".rds")))
    cat("RESULT:", label, "terms=", nrow(tab),
        "FDR<0.05=", sum(tab$p.adjust < 0.05, na.rm = TRUE), "\n")
    TRUE
  }, error = function(e) {
    message("FAILED ", label, ": ", conditionMessage(e))
    writeLines(conditionMessage(e), file.path(out, paste0(label, ".status.txt")))
    FALSE
  })
}
for (direction in c("increased", "decreased")) {
  genes <- if (direction == "increased") up else down
  run_one(paste0("GO_BP_ORA_", direction), function() clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500))
  run_one(paste0("Reactome_ORA_", direction), function() ReactomePA::enrichPathway(
    gene = genes, universe = universe, organism = "human", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500))
  run_one(paste0("KEGG_ORA_", direction), function() clusterProfiler::enrichKEGG(
    gene = genes, universe = universe, organism = "hsa", keyType = "ncbi-geneid",
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500))
}
run_one("GO_BP_GSEA", function() clusterProfiler::gseGO(
  geneList = rankvec, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
  ont = "BP", pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10,
  maxGSSize = 500, by = "fgsea", verbose = FALSE, seed = TRUE))
run_one("Reactome_GSEA", function() ReactomePA::gsePathway(
  geneList = rankvec, organism = "human", pvalueCutoff = 1,
  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
  by = "fgsea", verbose = FALSE, seed = TRUE))
run_one("KEGG_GSEA", function() clusterProfiler::gseKEGG(
  geneList = rankvec, organism = "hsa", keyType = "ncbi-geneid",
  pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10,
  maxGSSize = 500, by = "fgsea", verbose = FALSE, seed = TRUE))
writeLines(c(paste("R:", R.version.string),
             paste("Bioconductor:", as.character(BiocManager::version())),
             vapply(required, function(p) paste(p, as.character(packageVersion(p))), character(1))),
           file.path(out, "environment_versions.txt"))
cat("Phase 15E finished. Check .status.txt files and mapping_qa.tsv before interpretation.\n")
