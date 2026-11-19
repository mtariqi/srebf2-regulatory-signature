# Phase 16 shared helpers. Require explicit input provenance and fail closed.
options(stringsAsFactors=FALSE)
core_genes <- c('TM7SF2','CYP51A1','MSMO1','SQLE','FDFT1','HSD17B7',
                'LSS','NSDHL','DHCR7','EBP','SC5D','DHCR24')
read_tsv <- function(path) {
  if (!file.exists(path)) stop('Input absent: ', path)
  read.delim(path, sep='\t', check.names=FALSE, quote='', comment.char='', stringsAsFactors=FALSE)
}
choose_col <- function(d, choices, optional=FALSE) {
  hit <- intersect(choices, colnames(d))
  if (length(hit)) return(hit[1])
  if (optional) return(NULL)
  stop('Expected one of [', paste(choices, collapse=', '), ']; found [', paste(names(d), collapse=', '), ']')
}
map_symbols <- function(ids, supplied=NULL) {
  if (!is.null(supplied)) {
    s <- as.character(supplied)
    s[is.na(s) | !nzchar(s)] <- ids[is.na(s) | !nzchar(s)]
    return(s)
  }
  if (!requireNamespace('AnnotationDbi',quietly=TRUE) || !requireNamespace('org.Hs.eg.db',quietly=TRUE))
    stop('Ensembl IDs need Bioconductor AnnotationDbi + org.Hs.eg.db to map gene symbols')
  clean <- sub('\\.[0-9]+$', '', ids)
  isens <- grepl('^ENSG[0-9]+$',clean)
  symbols <- ids
  if (any(isens)) {
    mapping <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
      keys=unique(clean[isens]), column='SYMBOL', keytype='ENSEMBL', multiVals='first')
    mapped <- unname(mapping[clean[isens]])
    pos <- which(isens)
    symbols[pos[!is.na(mapped)]] <- mapped[!is.na(mapped)]
  }
  symbols
}
identify_gene_ids <- function(df) {
  col <- choose_col(df,c('gene_id','ensembl_gene_id','ENSEMBL','gene','Geneid','GeneID','id'),TRUE)
  if (is.null(col)) col <- colnames(df)[1]
  ids <- as.character(df[[col]])
  if (anyNA(ids) || any(!nzchar(ids))) stop('Missing gene IDs in ', col)
  list(column=col, ids=ids)
}
ensure_dir <- function(out) dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
