#!/usr/bin/env Rscript
# Real count matrix -> DESeq2 VST -> row Z scores; never treat raw counts as normalized.
suppressPackageStartupMessages({library(DESeq2); library(pheatmap); library(RColorBrewer)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop('Usage: Rscript phase16c_sterol_heatmap.R <raw_count_matrix.tsv> <sample_metadata.tsv> <output_prefix>')
source('scripts/phase16/phase16_utils.R')
counts <- read_tsv(args[1]); metadata <- read_tsv(args[2]); out <- args[3]; ensure_dir(out)
id <- identify_gene_ids(counts)
runcol <- choose_col(metadata,c('run','sample_id','sample','Sample'))
gcol <- choose_col(metadata,c('genotype','Genotype'))
tcol <- choose_col(metadata,c('treatment','condition','Treatment'))
run <- as.character(metadata[[runcol]])
if(length(run)!=12 || anyDuplicated(run)) stop('Expected 12 unique sample IDs in metadata')
if (!all(run %in% colnames(counts))) stop('Count matrix missing sample columns: ',paste(setdiff(run,colnames(counts)),collapse=', '))
mat <- as.matrix(counts[,run,drop=FALSE]);storage.mode(mat)<-'numeric'
if(any(!is.finite(mat))||any(mat<0) || any(abs(mat-round(mat))>1e-8)) stop('Expected finite nonnegative integer raw counts')
rownames(mat) <- id$ids
if(anyDuplicated(rownames(mat))) stop('Duplicate gene IDs')
mat <- round(mat)
# Genes with zero counts provide no information; include all 12 samples for normalization.
mat <- mat[rowSums(mat)>=10,,drop=FALSE]
metadata <- metadata[match(colnames(mat),metadata[[runcol]]),,drop=FALSE]
rownames(metadata) <- colnames(mat)
metadata[[gcol]] <- factor(metadata[[gcol]],levels=c('WT','KO'))
metadata[[tcol]] <- factor(metadata[[tcol]])
if(anyNA(metadata[[gcol]])|| anyNA(metadata[[tcol]])) stop('Unexpected genotype/treatment values')
coldata <- data.frame(genotype=metadata[[gcol]],treatment=metadata[[tcol]],row.names=rownames(metadata))
dds <- DESeq2::DESeqDataSetFromMatrix(countData=mat,colData=coldata,design=~ genotype + treatment)
dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::varianceStabilizingTransformation(dds,blind=TRUE)
vmat <- SummarizedExperiment::assay(vsd)
# Map Ensembl IDs to gene symbols, without silently choosing among duplicates.
sym <- map_symbols(rownames(vmat))
coreidx <- lapply(core_genes,function(g) which(sym==g))
if(any(lengths(coreidx)!=1L)) stop('Cannot map 12 sterol genes uniquely. Counts per gene: ',
      paste(paste(core_genes,lengths(coreidx),sep='='),collapse=', '))
idx <- unlist(coreidx,use.names=FALSE)
selected <- vmat[idx,,drop=FALSE]
rownames(selected) <- core_genes
zs <- t(scale(t(selected)))
if(any(!is.finite(zs))) stop('A sterol gene has zero variance; inspect VST values')
annotation <- coldata[,c('genotype','treatment'),drop=FALSE]
order <- order(annotation$genotype,annotation$treatment)
selected <- selected[,order,drop=FALSE];zs <- zs[,order,drop=FALSE];annotation<-annotation[order,,drop=FALSE]
cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(101)
makeplot <- function(filename,kind) {
 if(kind=='png') png(filename,width=2100,height=1300,res=220) else pdf(filename,width=10.5,height=6.5)
 pheatmap::pheatmap(zs,annotation_col=annotation,cluster_cols=FALSE,cluster_rows=FALSE,
   gaps_col=cumsum(table(interaction(annotation$genotype,annotation$treatment,drop=TRUE)))[-4],
   fontsize=10,fontsize_row=11,color=cols,breaks=seq(-2.5,2.5,length.out=102),
   main='HeLa SREBF2: 12-gene sterol core (VST, gene-wise Z score)',
   angle_col=45,border_color=NA)
 dev.off()
}
makeplot(paste0(out,'.png'),'png');makeplot(paste0(out,'.pdf'),'pdf')
write.table(data.frame(gene=rownames(selected),selected,check.names=FALSE),paste0(out,'_VST.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(gene=rownames(zs),zs,check.names=FALSE),paste0(out,'_row_Z.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(run=rownames(annotation),annotation),paste0(out,'_sample_order.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('PASS: HeLa sterol heatmap: ',nrow(zs),' genes, ',ncol(zs),' samples; all expression measured.\n')
