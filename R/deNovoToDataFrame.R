#' Creates data.frame with de novo HOMER results
#' 
#' Generates a data.frame containing the significant motifs found by HOMER.
#' @param results_homer Character string containing the path of the output folder that was generated 
#' with HOMER.
#' @param include_all TRUE (default) of FALSE for including also those motifs that have no ensembl_gene_id.
#' @return Data.frame containing the results and gene annotation for the de novo motifs.
#' @export
#' @examples
#' \dontrun{
#' df.tf <- deNovoToDataFrame("results_homer", include_all=FALSE)
#' }
deNovoToDataFrame <- function(results_homer, include_all=TRUE) {
  # Add "/" if not present in path
  if (!(grepl("*/$", results_homer))) results_homer=paste0(results_homer, "/")
  
  # Load data into data.frame
  res <- XML::readHTMLTable(paste0(results_homer, "homerResults.html"),
                            colClasses=c("character", "character", "numeric", "numeric",
                                         "character", "character", "character", "character", 
                                         "character"),
                            trim=TRUE, as.data.frame=TRUE, stringsAsFactors=FALSE)[[1]]
  fin <- grep("*", res$Rank, fixed=T)[1]-1 # Select number of significant motifs
  
  # Check if there are significant results
  if (fin==0) stop("There were no significant results. No data.frame generated.") 
  
  res <- res[1:fin,] # Select significant rows
  res$`% of Targets` <- as.numeric(gsub("%", "", res$`% of Targets`))
  res$`% of Background` <- as.numeric(gsub("%", "", res$`% of Background`))
  
  res$Motif.Name <- gsub("\\(.[[:digit:]]*.[[:digit:]]*.\\)More Information \\| Similar Motifs Found", 
                         "", 
                         res$`Best Match/Details`)
  
  if (include_all==TRUE) {
    res.all <- merge(res, dictHomerMotifs)
    not.inc <- res[!(res$Motif.Name %in% res.all$Motif.Name),]
    not.inc$Motif.Symbol <- unlist(lapply(strsplit(not.inc$Motif.Name, "/", fixed=TRUE), 
                                          function(x) x[[1]]))
    not.inc$external_gene_name <- NA
    not.inc$ensembl_gene_id <- NA
    not.inc <- not.inc[,c(10,1:9,11:13)]
    res.all <- rbind(res.all, not.inc)
  } else res.all <- merge(res, dictHomerMotifs)
  
  res.all$Motif.Score <- as.numeric(unlist(lapply(regmatches(res.all$`Best Match/Details`, 
                                    gregexpr("(?<=\\().*?(?=\\))", 
                                             res.all$`Best Match/Details`, perl=T)),
                                    function (x) x[[length(x)]])))
  res.all <- res.all[,c(2,14,4:7,1,11:13)]
  return(res.all)
}