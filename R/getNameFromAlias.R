#' Get gene symbol from alias.
#' 
#' Obtains a the gene symbol that corresponds to a gene alias.
#' @param alias Character vector containing the aliases for which you want to obtain the gene symbols.
#' @return Data.frame with the gene symbol and the alias.
#' @export
#' @examples 
#' \dontrun{
#' geneNames <- getNameFromAlias(alias)
#' }
getNameFromAlias <- function(alias) {
  queryGeneNames <- alias
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  # write your SQL query
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  # execute the query on the database
  aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
  # subset to get your results
  result <- aliasSymbol[toupper(aliasSymbol$alias_symbol) %in% toupper(as.character(queryGeneNames)),c(2,5)]
  return(result)
} 