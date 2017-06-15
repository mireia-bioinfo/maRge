#' Count significant results from de Novo motif analysis
#' 
#' Reads the html output from homer and counts how many motifs were classified as significant.
#' @param results_homer Either the path of the output folder that was generated containing HOMER reports or
#' the data.frame obtained using XML::readHTMLTable.
#' @return The number of significant motifs (int).
#' @export
#' @examples
#' \dontrun{
#' num <- countSignHOMER("results_homer")
#' }
countSignHOMER <- function(results_homer) {
  if (is(results_homer, "data.frame")) {
    tables <- results_homer
  } else {
  tables <- XML::readHTMLTable(paste0(results_homer, "homerResults.html"))[[1]]
  }
  fin <- grep("*", tables$Rank, fixed=T)[1]-1
  return(fin)
}