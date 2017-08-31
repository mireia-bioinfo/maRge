#' Creates data.frame with de novo HOMER results
#' 
#' Generates a data.frame containing the significant motifs found by HOMER.
#' @param results_homer Either the path of the output folder that was generated containing HOMER reports or
#' the data.frame obtained using XML::readHTMLTable.
#' @return The number of significant motifs (int).
#' @export
#' @examples
#' \dontrun{
#' num <- countSignHOMER("results_homer")
#' }
deNovoToDataFrame <- function(results_homer) {
  # Add "/" if not present in path
  if (!(grepl("*/$", results_homer))) results_homer=paste0(results_homer, "/")
  
  # Load data into data.frame
  res <- XML::readHTMLTable(paste0(results_homer, "homerResults.html"))[[1]]
  fin <- grep("*", res$Rank, fixed=T)[1]-1 # Select number of significant motifs
  
  # Check if there are significant results
  if (fin==0) stop("There were no significant results. No data.frame generated.") 
  
  res <- res[1:fin,] # Select significant rows
  res$`P-value` <- as.numeric(as.character(res$`P-value`))
  res$`log P-pvalue` <- as.numeric(as.character(res$`log P-pvalue`))
  res$`% of Targets` <- as.numeric(gsub("%", "", res$`% of Targets`))
  res$`% of Background` <- as.numeric(gsub("%", "", res$`% of Background`))
  
}