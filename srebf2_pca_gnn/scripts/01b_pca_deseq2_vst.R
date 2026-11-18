#!/usr/bin/env Rscript
# Recommended PCA: DESeq2 blind VST for sample-level QC, not inference.
suppressPackageStartupMessages({library(DESeq2)})
a <- commandArgs(trailingOnly=TRUE)
if(length(a)<3) stop('Usage: Rscript 01b_pca_deseq2_vst.R counts.tsv metadata.tsv outdir')
counts <- read.delim(a[1],check.names=FALSE)
md <- read.delim(a[2],check.names=FALSE)
if (!all(c('run','genotype','treatment') %in% colnames(md))) stop('Metadata missing columns')
idx <- match(md$run,colnames(counts))
if(anyNA(idx)) stop('Count matrix does not contain all sample IDs')
X <- as.matrix(counts[,idx,drop=FALSE])
if(anyNA(X)||any(X<0)||any(abs(X-round(X))>0)) stop('Expected nonnegative integer read counts')
storage.mode(X)<-'integer'
rownames(X)<-make.unique(as.character(counts[[1]]))
colnames(X)<-md$run
md$genotype<-factor(md$genotype,levels=c('WT','KO'))
md$treatment<-factor(md$treatment)
rownames(md)<-md$run
out<-a[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
for (excluded in c(NA,'SRR28966285')) {
  keep<- if(is.na(excluded)) rep(TRUE,nrow(md)) else md$run!=excluded
  suffix<-if(is.na(excluded)) 'all12' else 'exclude_SRR28966285'
  dds<-DESeqDataSetFromMatrix(X[,keep,drop=FALSE],md[keep,,drop=FALSE],design=~genotype*treatment)
  dds<-dds[rowSums(counts(dds)>=10)>=3,]
  vsd<-vst(dds,blind=TRUE)
  variability<-matrixStats::rowVars(assay(vsd))
  top<-head(order(variability,decreasing=TRUE),2000)
  p<-prcomp(t(assay(vsd)[top,,drop=FALSE]),center=TRUE,scale.=FALSE)
  scores<-data.frame(run=rownames(p$x),md[rownames(p$x),c('genotype','treatment'),drop=FALSE],PC1=p$x[,1],PC2=p$x[,2],row.names=NULL)
  write.table(scores,file=file.path(out,paste0('PCA_DESeq2_VST_',suffix,'_scores.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  vr<-p$sdev^2/sum(p$sdev^2)
  write.table(data.frame(PC=seq_along(vr),variance_fraction=vr),file=file.path(out,paste0('PCA_DESeq2_VST_',suffix,'_variance.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  for(ext in c('png','pdf')) {
    file<-file.path(out,paste0('PCA_DESeq2_VST_',suffix,'.',ext))
    if(ext=='png') png(file,width=2100,height=1600,res=250) else pdf(file,width=8.4,height=6.4)
    col<-ifelse(scores$genotype=='KO','#B24C42','#247090')
    pch<-ifelse(scores$treatment=='FBS',16,15)
    plot(scores$PC1,scores$PC2,col=col,pch=pch,cex=1.8,
         xlab=sprintf('PC1 (%.1f%%)',100*vr[1]),ylab=sprintf('PC2 (%.1f%%)',100*vr[2]),
         main=paste0('HeLa SREBF2 RNA-seq — DESeq2 VST (',suffix,')'))
    text(scores$PC1,scores$PC2,labels=substr(scores$run,nchar(scores$run)-3,nchar(scores$run)),pos=3,cex=.65)
    legend('topright',legend=c('KO','WT','FBS','depletion'),pch=c(16,16,16,15),col=c('#B24C42','#247090','#555555','#555555'),bty='n')
    dev.off()
  }
  message('PASS: ',suffix,'; VST PCA samples=',nrow(scores))
}
