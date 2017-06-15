#' Get annotation from TF name.
#' 
#' Obtains a summary description and the full name of the provided TF names.
#' @param namesTF Character string containing the names of the TF found by HOMER (it will remove "()" and "_".
#' @return Data.frame with the symbol of the TF/genes, their entrez id, ENSEMBL id, description and summary.
#' @export
#' @examples 
#' \dontrun{
#' anno <- getAnnoTF(namesTF=tf.names)
#' head(anno)
#'> external_gene_name TF.fullNames TF.symbol entrezgene ensembl_gene_id                       description                   summary
#'>               ATF3         Atf3      Atf3        467 ENSG00000162772 activating transcription factor 3   This gene encodes [...]
#' }
getAnnoTF <- function(namesTF) {
  namesTF <- unlist(sapply(sapply(namesTF, strsplit, "(", fixed=T), function(x) x[1]))
  namesTF <- unlist(sapply(sapply(namesTF, strsplit, "_", fixed=T), 
                           function(x) if (length(x) == 1) {x[1]} else {x[2]}))
  df <- data.frame("TF.fullNames"=names(namesTF),
                   "TF.symbol"=namesTF,
                   stringsAsFactors=F)
  ensembl = biomaRt::useEnsembl(biomart='ensembl', GRCh=37, dataset='hsapiens_gene_ensembl')
  dat.hs = biomaRt::getBM(attributes=c('external_gene_name', 'entrezgene'), filters = 'external_gene_name', 
                       values = df$TF.symbol, mart = ensembl)
  # Those that  were not found might be named after aliases
  alias <- df$TF.symbol[!(toupper(df$TF.symbol) %in% dat.hs$external_gene_name)]
  if (length(alias) > 0) {
    alias.df <- getNameFromAlias(as.character(alias))
    df$TF.symbol[toupper(df$TF.symbol) %in% toupper(alias.df$alias_symbol)] <- alias.df$symbol
    
    # Match alias with symbol
    for (i in 1:nrow(alias.df)) {
      df$TF.symbol[grep(alias.df$alias_symbol[i], df$TF.symbol, ignore.case=T)] <- alias.df$symbol[i]
    }
    
    # Run annotation again
    dat.hs = biomaRt::getBM(attributes=c('external_gene_name', 'entrezgene', 'ensembl_gene_id'), 
                            filters = 'external_gene_name', 
                            values = df$TF.symbol, mart = ensembl)
  }
  df$external_gene_name <- toupper(df$TF.symbol)
  df.anno <- merge(df, dat.hs, all.x=T)
  df.anno$description <- NA
  df.anno$summary <- NA
  
  for (i in df.anno$entrezgene) {
    r_search <- rentrez::entrez_search(db="gene", term=paste0("(", i, 
                                                              "[UID]) AND (Homo sapiens[ORGN])"))
    gene <- rentrez::entrez_fetch(db="gene", id=r_search$ids, rettype="xml", parsed=TRUE)
    gene_list <- XML::xmlToList(gene)
    df.anno$description[df.anno$entrezgene==i] <- gene_list$Entrezgene$Entrezgene_gene$`Gene-ref`$`Gene-ref_desc`
    if ("Entrezgene_summary" %in% names(gene_list$Entrezgene))
    df.anno$summary[df.anno$entrezgene==i] <- gene_list$Entrezgene$Entrezgene_summary
  }
  return(df.anno)
}