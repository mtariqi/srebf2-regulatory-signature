#!/usr/bin/env Rscript
# Phase 6F: exploratory overrepresentation of independently frozen SREBF2 reversal genes.
# Run from repository root. Requires msigdbr, unless MSIGDB_GMT points to a local GMT.
options(stringsAsFactors=FALSE)
root <- getwd()
d6 <- file.path(root,'results/validation/GSE147770/interaction_v1')
d5 <- file.path(root,'results/validation/GSE147770/direction_audit_v1')
out <- file.path(root,'results/validation/GSE147770/enrichment_v1')
if (dir.exists(out)) stop('Refusing to overwrite: ',out)
for (p in c(d6,d5)) if (!dir.exists(p)) stop('Missing directory: ',p)
read_tsv <- function(p) {if(!file.exists(p))stop('Missing file: ',p);read.delim(p,check.names=FALSE)}
wald_path <- file.path(d6,'interaction_Wald_v1.tsv')
audit_path <- file.path(d5,'SIGNATURE_DIRECTION_AUDIT_v1.tsv')
w <- read_tsv(wald_path); a <- read_tsv(audit_path)
req_w <- c('unversioned_id','padj'); req_a <- c('unversioned_id','direction','expected_reversal','interaction_FDR05','plus_log2FoldChange','minus_log2FoldChange')
if (!all(req_w %in% names(w)) || !all(req_a %in% names(a))) stop('Input schema changed; inspect column names')
if(anyDuplicated(w$unversioned_id)||anyDuplicated(a$unversioned_id))stop('Duplicate gene identifiers')
# Prespecified background: only genes with an evaluable interaction FDR, not the full genome.
universe <- unique(as.character(w$unversioned_id[!is.na(w$padj)]))
if(length(universe)!=14171L) stop('Expected 14171 evaluable interaction genes; got ',length(universe))
truth <- function(x) !is.na(x) & as.character(x) %in% c('TRUE','T','true','1')
a$selected <- truth(a$expected_reversal) & truth(a$interaction_FDR05)
sets <- list(positive=a$unversioned_id[a$selected & a$direction=='decreased_in_KO'],
             negative=a$unversioned_id[a$selected & a$direction=='increased_in_KO'],
             combined=a$unversioned_id[a$selected])
expected <- c(positive=131L,negative=50L,combined=181L)
for(n in names(sets)) {sets[[n]]<-intersect(unique(as.character(sets[[n]])),universe);if(length(sets[[n]])!=expected[[n]])stop('Unexpected selected gene count for ',n,': ',length(sets[[n]]))}
# MSigDB via msigdbr or explicitly supplied GMT. Use Ensembl IDs, never assume symbols are unique.
gmt <- Sys.getenv('MSIGDB_GMT','')
if(nzchar(gmt)) {
  if(!file.exists(gmt))stop('MSIGDB_GMT does not exist: ',gmt)
  lines <- readLines(gmt,warn=FALSE)
  entries <- strsplit(lines,'\t',fixed=FALSE)
  pathways <- lapply(entries,function(x) unique(x[-c(1,2)]))
  names(pathways) <- vapply(entries,`[`,character(1),1)
  pathway_collection <- rep('user_GMT',length(pathways))
  source_label <- paste0('GMT:',normalizePath(gmt))
  # GMT MUST contain unversioned human Ensembl gene IDs (ENSG...), not gene symbols.
  if(!any(grepl('^ENSG',unlist(pathways,use.names=FALSE))))stop('GMT genes must be unversioned Ensembl IDs (ENSG...), not symbols')
} else {
  if(!requireNamespace('msigdbr',quietly=TRUE))stop('Install msigdbr: Rscript -e "install.packages(\"msigdbr\", repos=\"https://cloud.r-project.org\")"')
  d <- msigdbr::msigdbr(species='Homo sapiens')
  required <- c('gs_name','ensembl_gene','gs_collection')
  if(!all(required %in% names(d)))stop('msigdbr columns unexpected: ',paste(names(d),collapse=', '))
  # Curated collections: Hallmark, Reactome and GO biological process.
  chosen <- d$gs_collection=='H' | (d$gs_collection=='C2' & grepl('REACTOME',d$gs_name,ignore.case=TRUE)) | (d$gs_collection=='C5' & grepl('GOBP|GO_BP',d$gs_name,ignore.case=TRUE))
  d <- d[!is.na(chosen)&chosen & !is.na(d$ensembl_gene),,drop=FALSE]
  d$ensembl_gene <- sub('\\.[0-9]+$','',as.character(d$ensembl_gene))
  pathway_collection <- tapply(as.character(d$gs_collection),as.character(d$gs_name),function(x)x[1])
  pathways <- split(as.character(d$ensembl_gene),as.character(d$gs_name))
  pathways <- lapply(pathways,unique)
  source_label <- paste0('msigdbr ',as.character(packageVersion('msigdbr')),'; MSigDB version ',if('db_version' %in% names(d))paste(unique(d$db_version),collapse=',')else 'not reported')
}
# Hypergeometric upper tail = one-sided overrepresentation, BH within each gene-set family and target component.
result <- list(); memberships <- list()
for (target_name in names(sets)) {
  target <- sets[[target_name]]; N <- length(universe); K <- length(target)
  for (i in seq_along(pathways)) {
    members <- intersect(unique(pathways[[i]]),universe)
    M <- length(members)
    if(M<10L || M>500L)next
    hits <- intersect(target,members); x <- length(hits)
    p <- phyper(x-1,M,N-M,K,lower.tail=FALSE)
    key <- paste(target_name,names(pathways)[i],sep='::')
    result[[key]] <- data.frame(component=target_name,collection=unname(pathway_collection[[names(pathways)[i]]]),pathway=names(pathways)[i],background_N=N,target_K=K,pathway_M=M,overlap=x,expected_overlap=K*M/N,fold_enrichment=if(M*K>0)x/(K*M/N) else NA_real_,pvalue=p,overlap_ensembl=paste(sort(hits),collapse=';'))
    if(x>0L) memberships[[key]] <- data.frame(component=target_name,pathway=names(pathways)[i],unversioned_id=hits)
  }
}
if(!length(result))stop('No gene sets with 10–500 genes overlap the background; check identifiers')
res <- do.call(rbind,result);rownames(res)<-NULL
res$padj_BH <- ave(res$pvalue,interaction(res$component,res$collection,drop=TRUE),FUN=function(z)p.adjust(z,method='BH'))
res <- res[order(res$component,res$padj_BH,res$pvalue),]
# Biological theme flags are descriptive only; the primary analysis is unbiased over all pathways.
res$theme_flag <- ifelse(grepl('CHOLESTEROL|STEROL|LIPID|FATTY_ACID|SREBP|ENDOPLASMIC_RETICULUM|STRESS',res$pathway,ignore.case=TRUE),'lipid_sterol_or_stress','other')
# Map gene IDs to available symbols if present in the Phase 6A mapping file.
map_path <- file.path(root,'results/validation/GSE147770/audit_v4/SIGNATURE_MAPPING_v4.tsv')
annot <- a[a$selected,c('unversioned_id','direction','plus_log2FoldChange','minus_log2FoldChange'),drop=FALSE]
annot$interaction_padj <- w$padj[match(annot$unversioned_id,w$unversioned_id)]
if(file.exists(map_path)) {
  mp<-read_tsv(map_path)
  symbol_col<-intersect(c('gene_symbol','Gene.symbol','symbol'),names(mp))
  if(length(symbol_col)) {
    mm<-mp[!duplicated(mp$unversioned_id),c('unversioned_id',symbol_col[1]),drop=FALSE]
    annot$gene_symbol<-mm[[symbol_col[1]]][match(annot$unversioned_id,mm$unversioned_id)]
  }
}
summary <- do.call(rbind,lapply(names(sets),function(n){z<-res[res$component==n,,drop=FALSE];data.frame(component=n,selected_genes=length(sets[[n]]),pathways_tested=nrow(z),pathways_FDR05=sum(z$padj_BH<.05),min_padj=if(nrow(z))min(z$padj_BH)else NA_real_)}))
inputs <- c(wald_path,audit_path,if(file.exists(map_path))map_path)
if(nzchar(gmt))inputs<-c(inputs,gmt)
dir.create(out,recursive=TRUE)
write.table(res,file.path(out,'PATHWAY_ENRICHMENT_ALL_v1.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
write.table(res[res$padj_BH<.05,,drop=FALSE],file.path(out,'PATHWAY_ENRICHMENT_FDR05_v1.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
write.table(res[res$theme_flag!='other',,drop=FALSE],file.path(out,'LIPID_STEROL_STRESS_PATHWAYS_v1.tsv'),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
write.table(annot,file.path(out,'REVERSAL_GENES_ANNOTATED_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(summary,file.path(out,'ENRICHMENT_SUMMARY_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(path=normalizePath(inputs),md5=unname(tools::md5sum(inputs))),file.path(out,'INPUT_CHECKSUMS_v1.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('Phase 6F v1: descriptive exploratory ORA','Source: ',source_label,'Background: 14171 genes with nonmissing Phase 6D Wald interaction padj','Targets: expected reversal AND interaction FDR<0.05; positive 131, negative 50, combined 181','Gene sets: 10 <= background-mapped genes <= 500','Test: one-sided hypergeometric','Correction: BH separately by component and collection','Caution: gene correlation, selection effects, and power can affect enrichment; no direct-binding inference',paste('R:',R.version.string)),file.path(out,'RUN_MANIFEST_v1.txt'))
cat('=== PHASE 6F ENRICHMENT SUMMARY ===\n');print(summary,row.names=FALSE)
cat('=== TOP 15 COMBINED PATHWAYS ===\n');print(head(res[res$component=='combined',c('collection','pathway','overlap','fold_enrichment','padj_BH')],15),row.names=FALSE)
cat('=== LIPID/STRESS PATHWAYS FDR<0.05 ===\n');print(head(res[res$theme_flag!='other' & res$padj_BH<.05,c('component','pathway','overlap','padj_BH')],20),row.names=FALSE)
cat('Saved: ',out,'\n')
