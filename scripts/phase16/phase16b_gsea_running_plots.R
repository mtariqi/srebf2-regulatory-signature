#!/usr/bin/env Rscript
# True GSEA running-score plots using exactly the stored Phase15E geneList and geneSets.
# Optional --eps0 reruns fgseaMultilevel over the same GO/Reactome pathway families.
suppressPackageStartupMessages({library(ggplot2);library(fgsea);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3 || length(args)>4) stop('Usage: Rscript phase16b_gsea_running_plots.R <GO_BP_GSEA.rds> <Reactome_GSEA.rds> <out_prefix> [--eps0]')
source('scripts/phase16/phase16_utils.R')
go <- readRDS(args[1]); re <- readRDS(args[2]); prefix <- args[3]; ensure_dir(prefix)
eps0 <- length(args)==4 && identical(args[4],'--eps0')
if(length(args)==4 && !eps0) stop('Unknown option: ',args[4])
get_details <- function(object,label) {
  if(!isS4(object) || !all(c('result','geneList','geneSets') %in% slotNames(object))) stop('Invalid clusterProfiler GSEA RDS: ', label)
  result <- as.data.frame(object@result)
  list(result=result,ranks=object@geneList,sets=object@geneSets,name=label)
}
objects <- list(GO=get_details(go,'GO'),Reactome=get_details(re,'Reactome'))
patterns <- list(
  list(db='GO',pattern='sterol biosynthetic process',label='GO sterol biosynthesis'),
  list(db='GO',pattern='cholesterol biosynthetic process',label='GO cholesterol biosynthesis'),
  list(db='Reactome',pattern='cholesterol biosynthesis',label='Reactome cholesterol biosynthesis')
)
results <- list(); plots <- list(); audit <- list(); eps0tables <- list()
for (db in names(objects)) {
 obj <- objects[[db]]
 ranks <- as.numeric(obj$ranks);names(ranks)<-names(obj$ranks)
 if(any(!is.finite(ranks)) || anyDuplicated(names(ranks))) stop('Invalid ranks: ',db)
 ranks <- sort(ranks,decreasing=TRUE)
 if(eps0) {
   cat('Running eps=0 across entire ',db,' GSEA family (may take time)...\n')
   set.seed(1516)
   tab <- fgsea::fgseaMultilevel(pathways=obj$sets,stats=ranks,minSize=10,maxSize=500,
                                eps=0,nPermSimple=10000)
   tab$leadingEdge <- vapply(tab$leadingEdge,function(x) paste(x,collapse='/'),character(1))
   write.table(as.data.frame(tab),paste0(prefix,'_',db,'_eps0_family.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
   eps0tables[[db]] <- tab
 }
 for (entry in patterns) {
   if (entry$db!=db) next
   desc <- as.character(obj$result$Description)
   matches <- which(tolower(trimws(desc)) == tolower(entry$pattern))
   if(length(matches)==0) matches <- which(grepl(entry$pattern,desc,ignore.case=TRUE,fixed=TRUE))
   # No near-match substitution; fail when term absent or ambiguous.
   if(length(matches)!=1) stop('Pathway match not unique for ',entry$label,'; matches: ',paste(desc[matches],collapse=' | '))
   hit <- obj$result[matches,,drop=FALSE]
   id <- as.character(hit$ID[1]); gene_set <- obj$sets[[id]]
   if (is.null(gene_set)) stop('No stored gene set for ',id)
   common <- intersect(as.character(gene_set),names(ranks))
   if (length(common)<10) stop('Too few mapped leading pathway genes: ',entry$label)
   p <- fgsea::plotEnrichment(common,ranks)+
      labs(title=entry$label,subtitle=paste0('Original NES ',format(round(hit$NES[1],3),nsmall=3),
        ' | original BH FDR ',format(signif(hit$p.adjust[1],3),scientific=TRUE)),
      x='Ranked genes (descending KO versus WT Wald statistic)',y='Running enrichment score')+
      theme_bw(base_size=10)+theme(plot.title=element_text(face='bold',size=11))
   plots[[length(plots)+1]] <- p
   audit[[length(audit)+1]] <- data.frame(database=db,pathway_id=id,description=desc[matches],
      original_NES=hit$NES[1],original_padj=hit$p.adjust[1],
      set_size=length(common),rank_length=length(ranks))
 }
}
fig <- patchwork::wrap_plots(plots,ncol=2)
ggsave(paste0(prefix,'.png'),fig,width=12.8,height=8.4,dpi=320,bg='white')
ggsave(paste0(prefix,'.pdf'),fig,width=12.8,height=8.4)
write.table(do.call(rbind,audit),paste0(prefix,'_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: plotted ',length(plots),' genuine GSEA running curves; eps=0 rerun = ',eps0,'\n')
cat('NOTE: eps=0 rerun is a numerical sensitivity only; report original and rerun FDR separately.\n')
