#!/usr/bin/env Rscript
# Phase 6D: independent HAP1 SREBF2 KO x serum interaction, reproducible audit.
suppressPackageStartupMessages(library(DESeq2))
options(stringsAsFactors=FALSE)
root <- getwd()
count_path <- file.path(root,'data/raw/validation/GSE147770/GSE147770_STAR_readCounts.txt.gz')
audit_path <- file.path(root,'results/validation/GSE147770/audit_v4/MATRIX_ROW_AUDIT_v4.tsv')
map_path <- file.path(root,'results/validation/GSE147770/audit_v4/SIGNATURE_MAPPING_v4.tsv')
out <- file.path(root,'results/validation/GSE147770/interaction_v1')
for(p in c(count_path,audit_path,map_path)) if(!file.exists(p)) stop('Missing input: ',p)
if(dir.exists(out)) stop('Refusing to overwrite existing output: ',out)
expected_sha <- '3f4a2307ede59793b77ac50c38fa519ece012c0a2128fe56ea4c0f44f0b7d598'
sha <- function(p) unname(tools::md5sum(p)) # for input audit use SHA256 from sha256sum below
sha256 <- function(p) {
  z <- suppressWarnings(system2('sha256sum', shQuote(p), stdout=TRUE, stderr=TRUE))
  if(length(z)!=1L || !grepl('^[a-f0-9]{64} ',z)) stop('sha256sum failed: ',p)
  sub(' .*','',z)
}
if(sha256(count_path)!=expected_sha) stop('Count file SHA256 mismatch')
m <- read.delim(gzfile(count_path),check.names=FALSE,quote='',na.strings=c('NA',''),colClasses='character')
a <- read.delim(audit_path,check.names=FALSE)
if(nrow(m)!=38344L || nrow(a)!=nrow(m)) stop('Matrix/audit row count mismatch')
if(!identical(as.character(m[['Ensembl.ID']]),as.character(a[['Ensembl.ID']]))) stop('Matrix/audit IDs differ')
if(!('eligible_complete_unique_gene' %in% names(a))) stop('Eligibility column absent')
eligible <- as.character(a[['eligible_complete_unique_gene']]) %in% c('TRUE','True','true','1')
if(sum(eligible)!=37919L) stop('Eligible gene count unexpected')
samples <- as.vector(outer(c('W','S'),paste0(rep(1:3,each=2),rep(c('plus','min'),3)),paste0))
# Explicit sample ordering, independent of file column order.
samples <- c(paste0('W',rep(1:3,each=2),rep(c('plus','min'),3)),paste0('S',rep(1:3,each=2),rep(c('plus','min'),3)))
if(!all(samples %in% names(m))) stop('Required sample columns absent')
ids <- sub('\\.[0-9]+$','',m[['Ensembl.ID']][eligible])
if(anyNA(ids) || anyDuplicated(ids)) stop('Eligible IDs are not unique')
values <- as.matrix(m[eligible,samples,drop=FALSE]); storage.mode(values) <- 'numeric'
if(anyNA(values) || any(!is.finite(values)) || any(values<0) || any(values!=floor(values))) stop('Invalid eligible counts')
if(any(values > .Machine$integer.max)) stop('Counts exceed R integer range')
counts <- matrix(as.integer(values),nrow=nrow(values),dimnames=list(ids,samples))
meta <- data.frame(sample=samples,genotype=ifelse(substr(samples,1,1)=='W','WT','SREBF2_KO'),serum=ifelse(grepl('plus$',samples),'plus_FBS','minus_FBS'),replicate=as.integer(sub('^[WS]([123]).*$','\\1',samples)))
meta$genotype <- factor(meta$genotype,levels=c('WT','SREBF2_KO'))
meta$serum <- factor(meta$serum,levels=c('plus_FBS','minus_FBS'))
rownames(meta) <- meta$sample
if(any(table(meta$genotype,meta$serum)!=3)) stop('Unexpected replication')
# A priori low-count filtering, independent of condition labels.
keep <- rowSums(counts>=10)>=3
if(sum(keep)<1000) stop('Unexpected low number of tested genes')
dds <- DESeqDataSetFromMatrix(countData=counts[keep,,drop=FALSE],colData=meta,design=~genotype+serum+genotype:serum)
dds <- DESeq(dds,quiet=TRUE)
coefs <- resultsNames(dds)
message('=== MODEL COEFFICIENTS ==='); print(coefs)
main <- grep('^genotype_SREBF2_KO_vs_WT$',coefs,value=TRUE)
interaction <- grep('^genotypeSREBF2_KO\\.serumminus_FBS$',coefs,value=TRUE)
if(length(main)!=1 || length(interaction)!=1) stop('Unexpected coefficient names: ',paste(coefs,collapse=', '))
extract <- function(x,contrast_name) {
  d <- as.data.frame(x)
  d$unversioned_id <- rownames(d)
  d$contrast <- contrast_name
  d[,c('unversioned_id','contrast','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
}
res_plus <- results(dds,name=main,alpha=.05)
res_minus <- results(dds,contrast=list(c(main,interaction)),alpha=.05)
res_interaction <- results(dds,name=interaction,alpha=.05)
res_lrt <- DESeq(dds,test='LRT',reduced=~genotype+serum,quiet=TRUE)
res_lrt_int <- results(res_lrt,alpha=.05)
res_list <- list(ko_vs_wt_plus_FBS=extract(res_plus,'ko_vs_wt_plus_FBS'),ko_vs_wt_minus_FBS=extract(res_minus,'ko_vs_wt_minus_FBS'),interaction_Wald=extract(res_interaction,'interaction_Wald'),interaction_LRT=extract(res_lrt_int,'interaction_LRT'))
# Signature enrichment among interaction genes: one-sided Fisher exact, within tested eligible genes.
map <- read.delim(map_path,check.names=FALSE)
required <- c('unversioned_id','mapping_status','direction')
if(!all(required %in% names(map))) stop('Mapping columns absent')
map <- map[map$mapping_status=='unique_complete',,drop=FALSE]
if(nrow(map)!=1065L || anyDuplicated(map$unversioned_id)) stop('Signature mapping unexpected')
tested <- rownames(dds)
inter <- res_interaction
sig <- !is.na(inter$padj) & inter$padj<.05
names(sig) <- rownames(inter)
universe <- names(sig)[!is.na(inter$padj)]
res_enrich <- list()
for(direction in c('decreased_in_KO','increased_in_KO','all')) {
  target <- if(direction=='all') map$unversioned_id else map$unversioned_id[map$direction==direction]
  target <- intersect(target,universe)
  bg <- setdiff(universe,target)
  a1 <- sum(sig[target]); b1 <- length(target)-a1
  c1 <- sum(sig[bg]); d1 <- length(bg)-c1
  mat <- matrix(c(a1,b1,c1,d1),nrow=2,byrow=TRUE)
  ft <- fisher.test(mat,alternative='greater')
  res_enrich[[direction]] <- data.frame(direction=direction,target_tested=length(target),background_tested=length(bg),target_significant=a1,background_significant=c1,odds_ratio=unname(ft$estimate),pvalue=ft$p.value)
}
enrich <- do.call(rbind,res_enrich)
enrich$padj_BH <- p.adjust(enrich$pvalue,method='BH')
# Report results with no claim of causal/direct binding.
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(meta,file.path(out,'SAMPLE_METADATA_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(n in names(res_list)) write.table(res_list[[n]],file.path(out,paste0(n,'_v1.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(enrich,file.path(out,'SIGNATURE_INTERACTION_ENRICHMENT_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(sample=samples,size_factor=sizeFactors(dds),library_size=colSums(counts),genes_detected=colSums(counts>0)),file.path(out,'SAMPLE_QC_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary <- data.frame(contrast=names(res_list),tested=nrow(dds),non_na_padj=vapply(res_list,function(x) sum(!is.na(x$padj)),integer(1)),significant_FDR_05=vapply(res_list,function(x) sum(x$padj<.05,na.rm=TRUE),integer(1)))
write.table(summary,file.path(out,'RESULT_SUMMARY_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
manifest <- c(paste('R:',R.version.string),paste('DESeq2:',as.character(packageVersion('DESeq2'))),paste('counts_sha256:',sha256(count_path)),paste('audit_sha256:',sha256(audit_path)),paste('mapping_sha256:',sha256(map_path)),paste('design: ~ genotype + serum + genotype:serum'),paste('reference genotype: WT; serum: plus_FBS'),paste('tested_genes:',nrow(dds)),paste('filter: counts >= 10 in >= 3 samples'),paste('wald_main:',main),paste('wald_interaction:',interaction),paste('interaction_lrt_reduced: ~ genotype + serum'),paste('no source files modified: TRUE'))
writeLines(manifest,file.path(out,'RUN_MANIFEST_v1.txt'))
message('=== PHASE 6D RESULT SUMMARY ===');print(summary,row.names=FALSE)
message('=== SIGNATURE INTERACTION ENRICHMENT ===');print(enrich,row.names=FALSE)
message('Saved: ',out)
