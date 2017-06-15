#' Obtain names from HOMER TF motif names.
#' 
#' Obtains human readable names (SYMBOLS) from HOMER de novo motif name matching.
#' @param results_homer Either the path of the output folder that was generated containing HOMER reports or
#' the data.frame obtained using XML::readHTMLTable.
#' @param num_sign Number of significant motifs that were found.
#' @return A vector with the names of the TF. 
#' @export
#' @examples 
#' \dontrun{
#' getSignNamesHOMER("results_homer", 2)
#' }
getSignNamesHOMER <- function(results_homer, num_sign) {
  if (is(results_homer, "data.frame")) {
    tables <- results_homer
    names <- as.character(tables[,8])
  } else if (is(results_homer, "character")) {
    if (length(results_homer) > 1) {
      names <- results_homer
    } else {
    tables <- XML::readHTMLTable(paste0(path_output, "/homerResults.html")[[1]])
    names <- as.character(tables[,8])
    }
  }

  names <- sapply(names, strsplit, "/", fixed=T)
  names <- unlist(sapply(names, function(x) x[1]))
  names(names) <- NULL
  return(names)
}