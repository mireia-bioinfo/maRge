#' Align TF around a SNP of interest
#' 
#' By using the PWM and the genomic locations of the motif found, it adds blanks to the PWM to make them be centered in the SNP of interest.
#' @param list Data.frame outputed from getCandidateTFs().
#' @param PFMatrixList Output of getMotifMatrixJASPAR() using matrix="PFM", with all the TF you want to query.
#' @importClassesFrom motifStack pcm pf
#' @export
alignTF_SNP <- function(list, PFMatrixList) {
  id <- list$ID
  PFM_id <- PFMatrixList[TFBSTools::ID(PFMatrixList) %in% as.character(id)]
  pcm_obj <- lapply(names(PFM_id), 
                    function(.ele, PFM_id){new("pcm", mat=TFBSTools::Matrix(PFM_id[[.ele]]), name=TFBSTools::name(PFM_id[[.ele]]),
                                               background=TFBSTools::bg(PFM_id[[.ele]]))}
                    ,PFM_id)
  names(pcm_obj) <- names(PFM_id)
  motifs<-lapply(pcm_obj, motifStack::pcm2pfm)
  
  lengths <- data.frame("ID"=names(motifs),
                     "mot.length"=sapply(motifs, function(x) ncol(x@mat)))
  list <- merge(list, lengths)

  list.pos <- list[list$strand=="+", c("ID", "external_gene_name", "relScore", "gen.start", "gen.end", "mot.length")]
  list.pos$offset.l <- list.pos$gen.start - min(list.pos$gen.start)
  list.pos$offset.r <- max(list.pos$gen.end) - list.pos$gen.end
  list.pos$check <- list.pos$offset.l + list.pos$mot.length + list.pos$offset.r
  
  list.pos$ID <- as.character(list.pos$ID)
  
  mot.al.pos <- list()
  for (i in 1:nrow(list.pos)) {
    pfm <-motifStack:::addBlank(motifs[[list.pos$ID[i]]], n=list.pos$offset.l[i], b=FALSE)
    pfm <- motifStack:::addBlank(pfm, n=list.pos$offset.r[i], b=TRUE)
    pfm@name <- paste0(pfm@name, "_", list.pos$relScore[i])
    mot.al.pos[[paste0(list.pos$ID[i], "_", list.pos$offset.l[i], "-", list.pos$offset.r[i])]] <- pfm
  }
  return(mot.al.pos)
}