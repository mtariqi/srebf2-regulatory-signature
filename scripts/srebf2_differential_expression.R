# Run from repository root with Rscript --vanilla.
prefix <- Sys.getenv("CONDA_PREFIX")
stopifnot(nzchar(prefix))
.libPaths(file.path(prefix, "lib", "R", "library"), include.site=FALSE)
suppressPackageStartupMessages(library(DESeq2))
set.seed(20261002)
out <- "results/de/discovery_v1"
if (dir.exists(out)) stop("Output directory exists; preserve it and choose a new output version.")
meta <- read.delim("metadata/GSE267018_analysis_samples.tsv", stringsAsFactors=FALSE)
tab <- read.delim("results/counts/GSE267018_raw_count_matrix.tsv", check.names=FALSE)
stopifnot(all(c("run", "sample_accession", "genotype", "treatment", "group") %in% names(meta)))
stopifnot(nrow(meta)==12L, !anyDuplicated(meta$run), !anyDuplicated(meta$sample_accession))
stopifnot(names(tab)[1]=="gene_id", !anyDuplicated(tab$gene_id), !anyNA(tab$gene_id))
stopifnot(setequal(names(tab)[-1], meta$run))
counts <- as.matrix(tab[,meta$run,drop=FALSE])
stopifnot(is.numeric(counts), all(is.finite(counts)), all(counts>=0), all(counts==floor(counts)), max(counts)<=.Machine$integer.max)
storage.mode(counts) <- "integer"
rownames(counts) <- tab$gene_id
stopifnot(all(meta$genotype %in% c("WT","KO")), all(meta$treatment %in% c("FBS","depletion")))
stopifnot(all(meta$group==paste(meta$genotype,meta$treatment,sep="_")))
meta$genotype <- factor(meta$genotype,levels=c("WT","KO"))
meta$treatment <- factor(meta$treatment,levels=c("FBS","depletion"))
rownames(meta) <- meta$run
stopifnot(all(table(meta$genotype,meta$treatment)==3L))
# Define a single expression filter using the full discovery cohort.
keep <- rowSums(counts>=10L)>=3L
stopifnot(sum(keep)>0L)
dir.create(out,recursive=TRUE)
write.table(data.frame(gene_id=rownames(counts),retained=keep),file.path(out,"GENE_FILTER.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_SESSION_INFO.txt"))
writeLines(c("Design: ~ genotype * treatment; reference levels WT and FBS.",
 "Filter: at least 10 counts in at least 3 of the original 12 samples; fixed in sensitivity fit.",
 "Primary: KO versus WT under depletion. Secondary: KO versus WT under FBS.",
 "Interaction: (KO-WT under depletion) minus (KO-WT under FBS), on log2 scale.",
 "Sensitivity: exclude SRR28966285; re-estimate normalization and dispersions.",
 "DESeq2 Wald tests; BH adjustment separately per contrast; alpha=0.05; default independent filtering and Cook's filtering.",
 "Unshrunk log2 fold changes; no fold-change cutoff; no automatic outlier exclusion or batch covariate.",
 "No signature selection or biological conclusions established by this script."),file.path(out,"METHODS.txt"))
summaries <- list()
all_results <- list()
for (fit in c("all12","without_SRR28966285")) {
  samples <- if(fit=="all12") meta$run else setdiff(meta$run,"SRR28966285")
  d <- DESeqDataSetFromMatrix(counts[keep,samples,drop=FALSE],meta[samples,,drop=FALSE],design=~genotype*treatment)
  d <- DESeq(d,parallel=FALSE)
  rn <- resultsNames(d)
  stopifnot(all(c("genotype_KO_vs_WT","genotypeKO.treatmentdepletion") %in% rn))
  dest <- file.path(out,fit); dir.create(dest)
  writeLines(rn,file.path(dest,"COEFFICIENTS.txt"))
  write.table(data.frame(run=samples,size_factor=sizeFactors(d)),file.path(dest,"SIZE_FACTORS.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  contrasts <- list(
    KO_vs_WT_depletion=results(d,contrast=list(c("genotype_KO_vs_WT","genotypeKO.treatmentdepletion")),alpha=0.05),
    KO_vs_WT_FBS=results(d,name="genotype_KO_vs_WT",alpha=0.05),
    interaction=results(d,name="genotypeKO.treatmentdepletion",alpha=0.05))
  all_results[[fit]] <- contrasts
  for (name in names(contrasts)) {
    r <- as.data.frame(contrasts[[name]])
    write.table(data.frame(gene_id=rownames(r),r),file.path(dest,paste0(name,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    sig <- !is.na(r$padj) & r$padj<0.05
    summaries[[length(summaries)+1L]] <- data.frame(fit=fit,contrast=name,genes=nrow(r),genes_with_pvalue=sum(!is.na(r$pvalue)),genes_with_padj=sum(!is.na(r$padj)),FDR05=sum(sig),up=sum(sig & r$log2FoldChange>0,na.rm=TRUE),down=sum(sig & r$log2FoldChange<0,na.rm=TRUE))
  }
  pdf(file.path(dest,"DISPERSIONS.pdf")); plotDispEsts(d); dev.off()
  saveRDS(d,file.path(dest,"dds.rds"))
}
summary <- do.call(rbind,summaries)
write.table(summary,file.path(out,"DE_SUMMARY.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for (name in names(all_results$all12)) {
  a <- as.data.frame(all_results$all12[[name]])
  b <- as.data.frame(all_results$without_SRR28966285[[name]])
  stopifnot(identical(rownames(a),rownames(b)))
  compare <- data.frame(gene_id=rownames(a),log2FC_all12=a$log2FoldChange,log2FC_sensitivity=b$log2FoldChange,padj_all12=a$padj,padj_sensitivity=b$padj,delta_log2FC=b$log2FoldChange-a$log2FoldChange)
  write.table(compare,file.path(out,paste0(name,"_SENSITIVITY.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
print(summary,row.names=FALSE)
cat("PASS: both fits and all three contrasts exported. Review diagnostics before interpretation.\n")
