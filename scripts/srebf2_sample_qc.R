# Run from repository root with Rscript --vanilla scripts/srebf2_sample_qc.R
prefix <- Sys.getenv('CONDA_PREFIX')
stopifnot(nzchar(prefix))
.libPaths(file.path(prefix, 'lib', 'R', 'library'), include.site=FALSE)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(pheatmap)})
out <- 'results/qc/cohort/expression'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- read.delim('results/counts/GSE267018_raw_count_matrix.tsv', check.names=FALSE)
stopifnot(!anyDuplicated(x$gene_id))
cts <- as.matrix(x[,-1]); rownames(cts) <- x$gene_id
stopifnot(is.numeric(cts), all(is.finite(cts)), all(cts >= 0), all(cts == floor(cts)))
storage.mode(cts) <- 'integer'
m <- read.csv('metadata/GSE267018_primary_runs.csv', check.names=FALSE)
stopifnot(!anyDuplicated(m$Run), setequal(m$Run, colnames(cts)))
m <- m[match(colnames(cts), m$Run), ]
stopifnot(all(m$genotype %in% c('Wild type','SREBP2-/-')),
          all(m$treatment %in% c('FBS','Cholesterol depletion')))
meta <- data.frame(run=m$Run, sample_accession=m[['Sample Name']],
 genotype=factor(ifelse(m$genotype=='Wild type','WT','KO'), levels=c('WT','KO')),
 treatment=factor(ifelse(m$treatment=='FBS','FBS','depletion'), levels=c('FBS','depletion')),
 row.names=m$Run)
meta$group <- factor(paste(meta$genotype, meta$treatment, sep='_'))
stopifnot(all(table(meta$group)==3), identical(rownames(meta),colnames(cts)))
write.table(meta, 'metadata/GSE267018_analysis_samples.tsv', sep='\t', quote=FALSE, row.names=FALSE)
# QC-only filter; record explicitly. No biological sample exclusions.
keep <- rowSums(cts >= 10) >= 3
stopifnot(sum(keep)>500)
dds <- DESeqDataSetFromMatrix(cts[keep,], meta, design=~genotype*treatment)
dds <- estimateSizeFactors(dds)
# Blind transformation for initial sample QC, not differential-expression testing.
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
a <- assay(vsd)
writeLines(c('Initial sample QC; all samples retained.',
 'Filter: at least 10 counts in at least 3 samples.',
 paste('Genes retained:',sum(keep)),
 'Transformation: DESeq2 varianceStabilizingTransformation, blind=TRUE.',
 'PCA: 500 genes with highest variance; center=TRUE, scale=FALSE.',
 'No differential-expression contrasts tested.'), file.path(out,'METHODS.txt'))
variance <- apply(a,1,var)
top <- order(variance,decreasing=TRUE)[seq_len(min(500,nrow(a)))]
pca <- prcomp(t(a[top,,drop=FALSE]),center=TRUE,scale.=FALSE)
percent <- 100*pca$sdev^2/sum(pca$sdev^2)
p <- cbind(meta, PC1=pca$x[,1], PC2=pca$x[,2])
write.table(p, file.path(out,'PCA_COORDINATES.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
plot <- ggplot(p,aes(PC1,PC2,color=genotype,shape=treatment))+geom_point(size=3.5)+
 geom_text(aes(label=sub('SRR28966','',run)),vjust=-0.8,size=3,show.legend=FALSE)+
 scale_color_manual(values=c(WT='#2878A5',KO='#D55E00'))+
 labs(x=sprintf('PC1 (%.1f%%)',percent[1]),y=sprintf('PC2 (%.1f%%)',percent[2]),
      caption='Labels show run accession suffixes. All 12 samples retained.')+
 theme_bw(base_size=12)+theme(legend.position='bottom')
ggsave(file.path(out,'PCA.png'),plot,width=9,height=7,dpi=300)
ggsave(file.path(out,'PCA.pdf'),plot,width=9,height=7)
correlation <- cor(a,method='pearson')
write.table(correlation,file.path(out,'SAMPLE_CORRELATIONS.tsv'),sep='\t',quote=FALSE,col.names=NA)
annotation <- meta[,c('genotype','treatment'),drop=FALSE]
pheatmap(correlation,annotation_col=annotation,annotation_row=annotation,
 filename=file.path(out,'SAMPLE_CORRELATIONS.png'),width=10,height=9,
 main='Pearson correlation of transformed counts',border_color=NA)
qc <- data.frame(run=colnames(cts),assigned_fragments=colSums(cts),
 nonzero_genes=colSums(cts>0),size_factor=sizeFactors(dds),group=meta$group)
write.table(qc,file.path(out,'LIBRARY_SIZE_QC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_SESSION_INFO.txt'))
cat('PASS: sample metadata, transformation, PCA and correlations generated.\n')
print(table(meta$group));print(qc)
