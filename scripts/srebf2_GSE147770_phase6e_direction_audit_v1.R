#!/usr/bin/env Rscript
# Phase 6E: descriptive directional audit of frozen signature against DESeq2 interaction.
# Input discovery signature, Phase 6A mapping, and Phase 6D DESeq2 tables; no refitting.
suppressPackageStartupMessages(library(utils))
root <- getwd()
base <- file.path(root, 'results/validation/GSE147770/interaction_v1')
out <- file.path(root, 'results/validation/GSE147770/direction_audit_v1')
if (dir.exists(out)) stop('Output directory exists; refusing overwrite: ', out)
if (!dir.exists(base)) stop('Missing Phase 6D directory: ', base)
files <- list.files(base, pattern='\\.(tsv|csv)$', full.names=TRUE, recursive=TRUE)
cat('=== AVAILABLE PHASE 6D TABLES ===\n'); print(basename(files))
if (!length(files)) stop('No Phase 6D TSV/CSV tables found')
readtab <- function(path) {
  if (grepl('\\.csv$', path, ignore.case=TRUE)) read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  else read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
}
# Identify DESeq2 tables by contents rather than guessing filenames.
classify <- function(f) {
  d <- tryCatch(readtab(f), error=function(e) NULL)
  if (is.null(d)) return(NULL)
  nm <- names(d)
  if (!all(c('log2FoldChange','padj') %in% nm)) return(NULL)
  list(path=f, columns=nm, n=nrow(d))
}
candidates <- Filter(Negate(is.null), lapply(files, classify))
cat('=== DESEQ2 TABLE CANDIDATES ===\n')
for (x in candidates) cat(basename(x$path), 'rows=',x$n,'columns=',paste(x$columns,collapse=','),'\n')
# Require explicit selection if naming is ambiguous; no silent interpretation.
choose <- function(pattern) {
  hits <- Filter(function(x) grepl(pattern,basename(x$path),ignore.case=TRUE), candidates)
  if (length(hits)!=1L) stop('Expected exactly one DESeq2 table matching ',pattern,'; found ',length(hits),'. Inspect list above and supply filenames to script.')
  hits[[1]]$path
}
plus_path <- choose('ko_vs_wt_plus_FBS')
minus_path <- choose('ko_vs_wt_minus_FBS')
interaction_path <- choose('interaction_Wald')
plus <- readtab(plus_path); minus <- readtab(minus_path); interaction <- readtab(interaction_path)
# DESeq2 results commonly store gene IDs as row names exported to unnamed first column.
get_id <- function(d, label) {
  options <- c('gene_id','unversioned_id','Ensembl.ID','ensembl_id','gene','X','...1','row.names','Row.names','')
  hits <- intersect(options,names(d))
  if (!length(hits)) stop(label, ': cannot locate gene identifier column. Found: ',paste(names(d),collapse=', '))
  v <- as.character(d[[hits[1]]]); v <- sub('\\.[0-9]+$','',v)
  if (anyNA(v) || anyDuplicated(v)) stop(label, ': missing or duplicate gene IDs')
  v
}
for (label in c('plus','minus','interaction')) {
  d <- get(label)
  d$unversioned_id <- get_id(d,label)
  assign(label,d)
}
map_path <- file.path(root,'results/validation/GSE147770/audit_v4/SIGNATURE_MAPPING_v4.tsv')
if (!file.exists(map_path)) stop('Missing frozen mapping: ',map_path)
map <- read.delim(map_path,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c('gene_id','unversioned_id','mapping_status','direction') %in% names(map)))
map <- map[map$mapping_status=='unique_complete',c('gene_id','unversioned_id','direction')]
if (anyDuplicated(map$unversioned_id)) stop('Duplicate signature IDs in mapping')
for (label in c('plus','minus','interaction')) {
  d <- get(label)
  d <- d[,c('unversioned_id','log2FoldChange','pvalue','padj'),drop=FALSE]
  names(d)[-1] <- paste0(label,'_',names(d)[-1])
  assign(label,d)
}
merged <- Reduce(function(x,y) merge(x,y,by='unversioned_id',all=TRUE),list(plus,minus,interaction))
merged <- merge(map,merged,by='unversioned_id',all.x=TRUE,sort=FALSE)
merged$sign_switch <- with(merged, !is.na(plus_log2FoldChange) & !is.na(minus_log2FoldChange) & plus_log2FoldChange * minus_log2FoldChange < 0)
merged$plus_negative <- !is.na(merged$plus_log2FoldChange) & merged$plus_log2FoldChange < 0
merged$minus_positive <- !is.na(merged$minus_log2FoldChange) & merged$minus_log2FoldChange > 0
merged$interaction_FDR05 <- !is.na(merged$interaction_padj) & merged$interaction_padj < 0.05
merged$expected_reversal <- ifelse(merged$direction=='decreased_in_KO', merged$plus_negative & merged$minus_positive, !merged$plus_negative & !merged$minus_positive & merged$sign_switch)
# More explicit: expected negative component is positive under +FBS, negative under -FBS.
merged$expected_reversal <- ifelse(merged$direction=='decreased_in_KO',
  !is.na(merged$plus_log2FoldChange) & !is.na(merged$minus_log2FoldChange) & merged$plus_log2FoldChange < 0 & merged$minus_log2FoldChange > 0,
  !is.na(merged$plus_log2FoldChange) & !is.na(merged$minus_log2FoldChange) & merged$plus_log2FoldChange > 0 & merged$minus_log2FoldChange < 0)
summary <- do.call(rbind,lapply(c('decreased_in_KO','increased_in_KO','all'),function(g){
  d <- if(g=='all') merged else merged[merged$direction==g,]
  data.frame(direction=g, mapped=nrow(d), both_LFC_available=sum(!is.na(d$plus_log2FoldChange)&!is.na(d$minus_log2FoldChange)),
  sign_switch=sum(d$sign_switch,na.rm=TRUE), expected_reversal=sum(d$expected_reversal,na.rm=TRUE),
  interaction_padj_available=sum(!is.na(d$interaction_padj)), interaction_FDR05=sum(d$interaction_FDR05),
  expected_reversal_and_interaction_FDR05=sum(d$expected_reversal & d$interaction_FDR05,na.rm=TRUE))
}))
# Wald vs LRT overlap: require LRT table, don't claim without it.
lrt_path <- choose('interaction_LRT')
lrt <- readtab(lrt_path); lrt$unversioned_id <- get_id(lrt,'LRT')
w <- readtab(interaction_path); w$unversioned_id <- get_id(w,'Wald')
overlap <- merge(w[,c('unversioned_id','padj')],lrt[,c('unversioned_id','padj')],by='unversioned_id',suffixes=c('_wald','_lrt'))
wald_sig <- !is.na(overlap$padj_wald)&overlap$padj_wald<0.05
lrt_sig <- !is.na(overlap$padj_lrt)&overlap$padj_lrt<0.05
concordance <- data.frame(wald_significant=sum(wald_sig),lrt_significant=sum(lrt_sig),both_significant=sum(wald_sig&lrt_sig),either_significant=sum(wald_sig|lrt_sig),jaccard=if(sum(wald_sig|lrt_sig)) sum(wald_sig&lrt_sig)/sum(wald_sig|lrt_sig) else NA_real_)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(merged,file.path(out,'SIGNATURE_DIRECTION_AUDIT_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(summary,file.path(out,'DIRECTION_SUMMARY_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(concordance,file.path(out,'WALD_LRT_OVERLAP_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
inputs <- c(map_path,plus_path,minus_path,interaction_path,lrt_path)
manifest <- data.frame(input=normalizePath(inputs),md5=unname(tools::md5sum(inputs)),stringsAsFactors=FALSE)
write.table(manifest,file.path(out,'INPUT_CHECKSUMS_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('=== DIRECTION SUMMARY ===\n'); print(summary,row.names=FALSE)
cat('=== WALD LRT OVERLAP ===\n'); print(concordance,row.names=FALSE)
cat('Saved:',out,'\n')
