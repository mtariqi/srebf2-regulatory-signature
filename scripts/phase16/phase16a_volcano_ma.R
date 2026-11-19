#!/usr/bin/env Rscript
# Validated all-12 DESeq2 primary cholesterol-depletion contrast; no invented values.
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript phase16a_volcano_ma.R <all12_KO_vs_WT_depletion.tsv> <output_prefix>')
source('scripts/phase16/phase16_utils.R')
d <- read_tsv(args[1]); output <- args[2]; ensure_dir(output)
id <- identify_gene_ids(d)
l2 <- choose_col(d,c('log2FoldChange','log2FC','log2FoldChange_shrunken'))
padj <- choose_col(d,c('padj','p.adjust','FDR'))
base <- choose_col(d,c('baseMean','mean_normalized_count','base_mean'))
sc <- choose_col(d,c('symbol','gene_symbol','SYMBOL','external_gene_name','hgnc_symbol'),TRUE)
genes <- map_symbols(id$ids, if(is.null(sc)) NULL else d[[sc]])
logfc <- as.numeric(d[[l2]]); q <- as.numeric(d[[padj]]); bm <- as.numeric(d[[base]])
if (anyDuplicated(id$ids)) stop('Duplicated gene IDs in DESeq2 output')
if (any(q<0 | q>1,na.rm=TRUE)) stop('Invalid FDR values')
if (any(bm<0,na.rm=TRUE)) stop('Negative baseMean')
positive <- q[is.finite(q) & q>0]
if(!length(positive)) stop('No positive adjusted p-values in DESeq2 results')
# Prevent zero-underflow p values from dominating the visual y axis.
qplot <- ifelse(!is.na(q) & q==0, min(positive)/10, q)
qplot <- pmax(qplot, .Machine$double.xmin)
classification <- ifelse(genes %in% core_genes,'Sterol core',
                     ifelse(!is.na(q)&q<0.05&is.finite(logfc)&abs(logfc)>=1,'FDR<0.05 and |LFC|>=1','Other'))
df <- data.frame(gene_id=id$ids,gene_symbol=genes,log2FoldChange=logfc,padj=q,
                 baseMean=bm,neglog10_FDR=-log10(qplot),class=classification)
df <- df[is.finite(df$log2FoldChange)&is.finite(df$baseMean)&!is.na(df$padj),]
df$class <- factor(df$class,levels=c('Other','FDR<0.05 and |LFC|>=1','Sterol core'))
if(!all(core_genes %in% df$gene_symbol)) warning('Sterol genes missing: ', paste(setdiff(core_genes,df$gene_symbol),collapse=', '))
labels <- df[df$gene_symbol %in% core_genes,]
cols <- c('Other'='grey75','FDR<0.05 and |LFC|>=1'='#BA5147','Sterol core'='#155E8B')
baseplot <- list(scale_color_manual(values=cols,drop=FALSE),theme_bw(base_size=11),
  theme(legend.position='bottom',plot.title=element_text(face='bold')))
p1 <- ggplot(df,aes(log2FoldChange,neglog10_FDR,color=class))+geom_point(alpha=.65,size=.9)+
  geom_vline(xintercept=c(-1,1),linetype=3,color='grey45')+
  geom_hline(yintercept=-log10(.05),linetype=3,color='grey45')+
  geom_text_repel(data=labels,aes(label=gene_symbol),color='#17324C',size=2.8,max.overlaps=Inf,seed=16)+
  labs(title='A  Volcano: KO versus WT',x='log2 fold change (KO / WT)',y='−log10(BH-adjusted p)',color=NULL)+baseplot
p2 <- ggplot(df,aes(log10(baseMean+1),log2FoldChange,color=class))+geom_point(alpha=.65,size=.9)+
  geom_hline(yintercept=0,color='grey40')+
  geom_hline(yintercept=c(-1,1),linetype=3,color='grey45')+
  geom_text_repel(data=labels,aes(label=gene_symbol),color='#17324C',size=2.8,max.overlaps=Inf,seed=16)+
  labs(title='B  MA: KO versus WT',x='log10(baseMean + 1)',y='log2 fold change (KO / WT)',color=NULL)+baseplot
fig <- p1+p2+plot_layout(guides='collect') & theme(legend.position='bottom')
ggsave(paste0(output,'.png'),fig,width=13,height=5.5,dpi=320,bg='white')
ggsave(paste0(output,'.pdf'),fig,width=13,height=5.5)
write.table(df,paste0(output,'_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS volcano/MA: ',nrow(df),' measured genes; sterol genes labeled ', nrow(labels),'\n')
