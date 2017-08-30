#' Get JASPAR2016 Motif Matrix
#' 
#' Obtains a PWMatrixList object containing all the PWM that match the query generated to JASPAR2016 database.
#' @param species Integer containing the taxa id for the sepecie you want to query. Default 9606 (Homo sapiens).
#' @param matrixtype Character containing the type of motif matrix you want to retrieve from JASPAR database. One of "PWM" (default), "ICM", "PFM". 
#' @return Object of type "PFMatrixList" that can be used as input for getCandidateTFs.
#' @export
#' @importFrom RSQLite dbGetQuery dbConnect dbDisconnect
#' @examples
#' \notrun{
#' getMotifMatrixJASPAR(species=9606, matrixtype="PWM")
#' }
getMotifMatrixJASPAR <- function(species=9606, matrixtype="PWM") {
  opts <- list()
  opts[["species"]] <- species
  opts[["matrixtype"]] <- matrixtype
  opts[["all_versions"]] <- TRUE
  PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts)
  return(PFMatrixList)
}