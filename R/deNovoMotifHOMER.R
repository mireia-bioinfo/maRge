#' HOMER de novo motif analysis 
#'
#' Finds instances of de Novo motifs in your regions of interest using HOMER suit.
#' @param bed Bed file containing the regions you want to analyze.
#' @param genome Name of the reference genome to use. "hg19" is set as default.
#' @param path_output Path of the output folder that will be generated containing HOMER reports.
#' @param size Corresponds to the flag "-size" in fingMotifsGenome.pl. Default is "given".
#' @param bits TRUE (default) or FALSE for scaling the logo.
#' @param other_param Character string containing all other parameters that should be added to the command line.
#' @param cores Integer indicating the number of cores to use for the analysis.
#' @param path_homer Path for homer location in your computer. It automatically looks in the aliases.
#' @param run TRUE (default) or FALSE for wether to run the command. If FALSE, will only print the line.
#' @return Folder "path_ouput" with the html resulting from HOMER de novo motif analysis (no R object).
#' @export
#' @examples 
#' \dontrun{
#' deNovoMotifHOMER(bed=positions, path_output=test_folder)
#' }
deNovoMotifHOMER <- function(bed, genome="hg19", path_output, cores=6, path_homer="default",
                             size="given", bits=TRUE, other_param="", run=TRUE) {
  ## Set HOMER path
  if (path_homer=="default") {
    path_homer=gsub("bin/homer", "", Sys.which("homer"), fixed=T)
    if (path_homer == "/") {
      stop("Could not find HOMER installation. Please provide path using argument 'path_homer'")
    }
  }
  
  ## Create command line for running homer
  homer <- paste(paste0(path_homer, "bin/findMotifsGenome.pl"), 
                 bed, genome, path_output,
                 paste("-p", cores),
                 "-size", size,
                 other_param
                 )
  if (bits==TRUE) homer <- paste0(homer, " -bits")
  
  print(homer) # Echo de line of code as message
  if (run==TRUE) system(homer) # Run line of code
}