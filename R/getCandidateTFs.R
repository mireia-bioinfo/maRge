#' Obtain candidate TF that might be binding a specified region
#' 
#' Using an input DNA sequence or genomic coordinates, it returns a list of candidate TF that match motif sequences found in your region.
#' @param dat Object of type GRanges or data.frame for genomic region (it will only use the first row of these objects), or character or DNAstring for type sequence.
#' @param PFMatrixList Object of type "PFMatrixList" containing the subset of motif matrices against whic compare our sequence.
#' @param region Character string with the type of dat you are using as input: "region" (default) for using genomic coordinates as input or "seq" for using a nucleotide sequence as input.
#' @param seqname Character string with the name of the region.
#' @param minscore Minimum percentage for matching motifs. Default: "80%".
#' @param strand String containing the strand in which to perform the analysis. "*" both strands (default), "+" positive strand and "-" negative strand.
#' @param start Integer containing the starting coordinates of region, if not provided using a data.frame or GRanges object.
#' @param end Character string containing chromosome name, if not provided using a data.frame or GRanges object.
#' @return Data.frame containing all the TF of interest that can be found in the region.
#' @export
#' @examples 
#' \notrun{
#' 
#' }
getCandidateTFs <- function(dat, PFMatrixList, type="region", seqname, 
                            minscore="80%", strand="*", start=NA, chr=NA) {
  ## Check that input has correct format.
  if (type=="region") {
    if (class(dat) == "GRanges") {
      seq <- Biostrings::DNAString(toString(XVector::subseq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19[[as.character(GenomicRanges::seqnames(dat))]], 
                                                            start=data.frame(GenomicRanges::ranges(dat))[1,1], 
                                                            end=data.frame(GenomicRanges::ranges(dat))[1,2])))
    }
    else if (class(dat)=="data.frame") {
      seq <- Biostrings::DNAString(toString(XVector::subseq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19[[as.character(dat[1,1])]], 
                                                            start=dat[1,2], 
                                                            end=dat[1,3])))
    }
    else {
      stop("Invalid dat object for type='region'. Should be a GRanges object or a data.frame (only first row will be used)", 
           call. = TRUE, domain = NULL)
    }
  }
  if (type=="seq") {
    if (class(dat)=="character") {
      seq <- Biostrings::DNAString(dat)
    }
    else if (class(dat)=="DNAString") {
      seq <- dat
    }
    else {
      stop("Invalid dat object for type='seq'. Should be a character or a DNAstring.", 
           call. = TRUE, domain = NULL)
    }
  }
  # else {
  #   stop("Invalid type specification. Should be 'region' or 'seq'.", 
  #        call. = TRUE, domain = NULL)
  # }
  
  sitesetList <- TFBSTools::searchSeq(PFMatrixList, seq, seqname=seqname,
                           min.score="80%", strand="*")
  list <- data.frame(sitesetList)[,-c(2,3)]
  
  if (type=="seq" & is.na(start) & is.na(chr)) {
    return(list)
  }
  if (type=="region") {
    if (class(dat) == "GRanges") {
      list$gen.chr <- as.character(GenomicRanges::seqnames(dat))
      list$gen.start <- data.frame(GenomicRanges::ranges(dat))[1,1] + list$start
      list$gen.end <- data.frame(GenomicRanges::ranges(dat))[1,1] + list$end
      return(list)
    }
    else if (class(dat)=="data.frame") {
      list$gen.chr <- as.character(dat[1,1])
      list$gen.start <- dat[1,2] + list$start
      list$gen.end <- dat[1,2] + list$end
      return(list)
    }
  }
  if (type=="seq" & !is.na(start) & !is.na(chr)){
    list$gen.chr <- as.character(chr)
    list$gen.start <- start + list$start
    list$gen.end <- start + list$end
    return(list)
  }
}