#' Join significant de novo motifs found by HOMER
#' 
#' Creates a unique .motif file containing all the significant motifs that were found
#' by HOMER de novo motif analysis.
#' @param path_output Path of the output folder that was generated containing HOMER reports.
#' @param num_sign Number of significant motifs that were found.
#' @return A file "path_output.motif" with all the motifs found significant in the regions.
#' @export
#' @examples 
#' \dontrun{
#' catSignMotifsHOMER("results_homer")
#' }
catSignMotifsHOMER <- function(path_output, num_sign) {
  loc <- paste0(path_output, "/homerResults/motif", 1:num_sign, ".motif")
  cmd <- paste("cat", paste(loc, collapse=" "), ">", 
               paste0(substr(path_output, 1, nchar(path_output)-1), ".motif"))
  system(cmd)
  return(invisible(paste0(substr(path_output, 1, nchar(path_output)-1), ".motif")))
}