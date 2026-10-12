#!/usr/bin/env Rscript
# Phase 6G: mitochondrial leading genes and enrichment robustness.
# Requires msigdbr; uses only frozen Phase 6D/6E/6F outputs.
options(stringsAsFactors=FALSE)
root <- getwd()
d6 <- file.path(root,'results/validation/GSE147770/interaction_v1')
d5 <- file.path(root,'results/validation/GSE147770/direction_audit_v1')
df <- file.path(root,'results/validation/GSE147770/enrichment_v1')
out <- file.path(root,'results/validation/GSE147770/mitochondrial_robustness_v1')
if(dir.exists(out)) stop('Output already exists (refusing overwrite): ',out)
paths <- c(wald=file.path(d6,'interaction_Wald_v1.tsv'),plus=file.path(d6,'ko_vs_wt_plus_FBS_v1.tsv'),minus=file.path(d6,'ko_vs_wt_minus_FBS_v1.tsv'),audit=file.path(d5,'SIGNATURE_DIRECTION_AUDIT_v1.tsv'),enrichment=file.path(df,'PATHWAY_ENRICHMENT_ALL_v1.tsv'))
for(p in paths) if(!file.exists(p))stop('Missing required input: ',p)
if(!requireNamespace('msigdbr',quietly=TRUE))stop('Install msigdbr first: Rscript -e "install.packages(\"msigdbr\",repos=\"https://cloud.r-project.org\")"')
rd <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
w <- rd(paths['wald']); plus <- rd(paths['plus']); minus <- rd(paths['minus']); a <- rd(paths['audit']); enr <- rd(paths['enrichment'])
for(z in list(w,plus,minus))if(!all(c('unversioned_id','log2FoldChange','padj')%in%names(z)))stop('DESeq2 input schema mismatch')
if(!all(c('unversioned_id','direction','expected_reversal','interaction_FDR05')%in%names(a)))stop('Audit schema mismatch')
if(!all(c('component','pathway','overlap','padj_BH')%in%names(enr)))stop('Phase 6F enrichment schema mismatch')
for(z in list(w,plus,minus,a))if(anyDuplicated(z$unversioned_id))stop('Duplicate gene identifiers in inputs')
truth <- function(z) !is.na(z) & as.character(z)%in%c('TRUE','true','T','1')
selected <- truth(a$expected_reversal)&truth(a$interaction_FDR05)
universe <- unique(as.character(w$unversioned_id[!is.na(w$padj)]))
if(length(universe)!=14171L)stop('Expected 14171 tested genes; found ',length(universe))
sets <- list(positive=intersect(a$unversioned_id[selected & a$direction=='decreased_in_KO'],universe),negative=intersect(a$unversioned_id[selected & a$direction=='increased_in_KO'],universe),combined=intersect(a$unversioned_id[selected],universe),all_interaction_signature=intersect(a$unversioned_id[truth(a$interaction_FDR05)],universe))
expected <- c(positive=131L,negative=50L,combined=181L,all_interaction_signature=359L)
for(n in names(sets))if(length(sets[[n]])!=expected[[n]])stop('Unexpected count for ',n,': ',length(sets[[n]]))
# Database version/identity must be recorded; Ensembl IDs are required.
db <- msigdbr::msigdbr(species='Homo sapiens')
if(!all(c('gs_name','ensembl_gene','gs_collection')%in%names(db)))stop('msigdbr schema mismatch')
db <- db[!is.na(db$ensembl_gene),,drop=FALSE]
db$ensembl_gene <- sub('\\.[0-9]+$','',as.character(db$ensembl_gene))
get_pathway <- function(name) unique(db$ensembl_gene[db$gs_name==name])
hallmark <- 'HALLMARK_OXIDATIVE_PHOSPHORYLATION'
hall <- get_pathway(hallmark)
if(!length(hall))stop('Hallmark OXPHOS absent in msigdbr; inspect database version')
# Verify database membership reproduces Phase 6F overlap exactly; avoid silent version drift.
leading <- sort(intersect(sets$combined,hall))
prior <- enr[enr$component=='combined' & enr$pathway==hallmark,,drop=FALSE]
if(nrow(prior)!=1L || length(leading)!=as.integer(prior$overlap[1]))stop('MSigDB membership/version differs from Phase 6F; do not compare without freezing database. Current overlap ',length(leading),' versus prior ',paste(prior$overlap,collapse=','))
if(length(leading)!=13L)stop('Expected 13 Hallmark leading genes; found ',length(leading))
# Phase 6D contrast coefficients; gene-level log2FC, not expression or protein activity.
get_fc <- function(tab,col) {x <- tab[,c('unversioned_id','log2FoldChange','padj')];names(x)[2:3]<-paste0(col,c('_log2FC','_padj'));x}
annot <- Reduce(function(x,y)merge(x,y,by='unversioned_id',all.x=TRUE,sort=FALSE),list(data.frame(unversioned_id=leading),get_fc(plus,'plus'),get_fc(minus,'minus'),get_fc(w,'interaction')))
annot <- merge(annot,a[,intersect(c('unversioned_id','direction','gene_id'),names(a)),drop=FALSE],by='unversioned_id',all.x=TRUE,sort=FALSE)
annot <- annot[match(leading,annot$unversioned_id),]
annot$effect_difference_minus_vs_plus <- annot$minus_log2FC-annot$plus_log2FC
annot$sign_switch <- annot$plus_log2FC*annot$minus_log2FC<0
# Gene symbol annotation is optional and descriptive; IDs remain authoritative.
if('Gene.symbol' %in% names(db)) { } # msigdbr uses gene_symbol, not Gene.symbol
if('gene_symbol'%in%names(db)){
  sy <- unique(db[db$gs_name==hallmark,c('ensembl_gene','gene_symbol')]);sy <- sy[!duplicated(sy$ensembl_gene),]
  annot$gene_symbol <- sy$gene_symbol[match(annot$unversioned_id,sy$ensembl_gene)]
}
# Per-gene OXPHOS pathway annotations and broad genomic origin check (ENSG IDs alone do not encode genome).
# Ensembl IDs do not determine mitochondrial-encoded status. Leave unknown without authoritative chromosome annotation.
annot$genomic_origin <- 'not_determined_from_available_inputs'
# Primary + alternative universe: signature genes with evaluable interaction padj.
signature_universe <- intersect(a$unversioned_id,universe)
if(length(signature_universe)!=807L)stop('Expected 807 evaluable signature genes, got ',length(signature_universe))
# Fisher/hypergeometric ORA; pathway membership and selected genes both restricted to each universe.
ora <- function(target,background,pathway,component,universe_label){
  target<-intersect(unique(target),background);members<-intersect(unique(pathway),background)
  N<-length(background);K<-length(target);M<-length(members);x<-length(intersect(target,members))
  data.frame(component=component,universe=universe_label,background_N=N,target_K=K,pathway_M=M,overlap=x,expected_overlap=K*M/N,fold_enrichment=if(K*M>0)x*N/(K*M)else NA_real_,pvalue=phyper(x-1,M,N-M,K,lower.tail=FALSE),overlap_ids=paste(sort(intersect(target,members)),collapse=';'))
}
# Enrichment robustness is evaluated for the Hallmark term only, so p-values here are nominal.
rob <- do.call(rbind,lapply(names(sets),function(n)do.call(rbind,lapply(list(all_evaluable=universe,signature_only=signature_universe),function(u){ora(sets[[n]],u,hall,n,if(identical(u,universe))'all_evaluable'else'signature_only')}))))
# Multiple-testing sensitivity across same 5647 pathways using Phase 6F results; show selected pathways and correction scope.
ox_related <- enr[enr$component=='combined' & grepl('OXIDATIVE_PHOSPHORYLATION|RESPIRATORY_ELECTRON_TRANSPORT|ELECTRON_TRANSPORT_CHAIN|ATP_SYNTHESIS|MITOCHONDRIAL_TRANSLATION',enr$pathway,ignore.case=TRUE),c('collection','pathway','overlap','fold_enrichment','pvalue','padj_BH')]
ox_related <- ox_related[order(ox_related$padj_BH),]
# A second full-universe comparison: same ORA for the 359-gene interaction-positive signature subset.
# This is a sensitivity analysis, not an independent validation.
# Make publication-friendly fold-change figure and compact heatmap without pretending log2FC is expression.
dir.create(out,recursive=TRUE)
write.table(annot,file.path(out,'OXPHOS_13_GENES_v1.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(rob,file.path(out,'OXPHOS_BACKGROUND_SENSITIVITY_v1.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
write.table(ox_related,file.path(out,'OXPHOS_RELATED_PATHWAYS_v1.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')
pdf(file.path(out,'OXPHOS_13_LFC_HEATMAP_v1.pdf'),width=7.8,height=7.8)
mat <- as.matrix(annot[,c('plus_log2FC','minus_log2FC')]);rownames(mat)<-if('gene_symbol'%in%names(annot))ifelse(is.na(annot$gene_symbol)|annot$gene_symbol=='',annot$unversioned_id,annot$gene_symbol)else annot$unversioned_id
colnames(mat)<-c('KO vs WT (+FBS)','KO vs WT (-FBS)')
lim<-max(abs(mat),na.rm=TRUE); if(!is.finite(lim)||lim==0)lim<-1
pal<-colorRampPalette(c('#2166ac','white','#b2182b'))(101)
par(mar=c(6,12,3,3)); image(x=1:2,y=1:nrow(mat),z=t(mat[nrow(mat):1,,drop=FALSE]),col=pal,zlim=c(-lim,lim),xaxt='n',yaxt='n',xlab='',ylab='',main='HAP1 SREBF2 KO effect: log2 fold changes')
axis(1,at=1:2,labels=colnames(mat),las=2,cex.axis=.85);axis(2,at=1:nrow(mat),labels=rev(rownames(mat)),las=2,cex.axis=.8);box();dev.off()
pdf(file.path(out,'OXPHOS_13_CONDITION_LFC_v1.pdf'),width=8.5,height=6.8)
ord <- order(annot$effect_difference_minus_vs_plus);plot(NA,xlim=c(0.5,2.5),ylim=range(c(annot$plus_log2FC,annot$minus_log2FC),finite=TRUE),xaxt='n',xlab='',ylab='KO vs WT log2 fold change',main='Hallmark OXPHOS reversal genes')
axis(1,at=1:2,labels=c('+FBS','-FBS'));abline(h=0,lty=2,col='gray50')
for(i in ord)lines(c(1,2),c(annot$plus_log2FC[i],annot$minus_log2FC[i]),col=adjustcolor('steelblue',alpha.f=.65),lwd=1.5)
points(rep(1,nrow(annot)),annot$plus_log2FC,pch=16,cex=.65);points(rep(2,nrow(annot)),annot$minus_log2FC,pch=16,cex=.65)
dev.off()
inputs<-unname(paths)
write.table(data.frame(path=normalizePath(inputs),md5=unname(tools::md5sum(inputs))),file.path(out,'INPUT_CHECKSUMS_v1.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c('Phase 6G v1; no DESeq2 refitting','Hallmark pathway: HALLMARK_OXIDATIVE_PHOSPHORYLATION','Primary Phase 6F background: 14171 genes with interaction padj','Sensitivity background: 807 evaluable frozen signature genes','Sensitivity selected target: 359 interaction-significant signature genes','One-sided hypergeometric ORA for Hallmark only; nominal p-values; no multi-pathway adjusted p from this sensitivity analysis',paste('msigdbr version:',as.character(packageVersion('msigdbr'))),paste('MSigDB version:',if('db_version'%in%names(db))paste(unique(db$db_version),collapse=',')else'unknown'),paste('R version:',R.version.string),'Genomic origin unknown without chromosome annotation; do not infer from ENSG ID','Figures display gene-level KO vs WT log2FC, not expression or protein function','Phase 6F enrichment adjusted p-values computed within component x collection; see prior manifest'),file.path(out,'RUN_MANIFEST_v1.txt'))
cat('=== PHASE 6G LEADING OXPHOS GENES ===\n');print(annot[,intersect(c('gene_symbol','unversioned_id','direction','plus_log2FC','minus_log2FC','interaction_padj'),names(annot))],row.names=FALSE)
cat('=== PHASE 6G OXPHOS BACKGROUND SENSITIVITY ===\n');print(rob[,c('component','universe','background_N','target_K','pathway_M','overlap','fold_enrichment','pvalue')],row.names=FALSE)
cat('=== PHASE 6G TOP RELATED PATHWAYS ===\n');print(head(ox_related,12),row.names=FALSE)
cat('Saved: ',out,'\n')
