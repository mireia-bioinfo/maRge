#' Draw list of motif logos in stack
#' 
#' Using a list of pfm objects, it draws a list of motifs and highlights the position where the SNP can be found.
#' @param list List of pfm objects. Also output from alignTF_SNP().
#' @param length Total length of the motifs that are being plotted.
#' @param posSNP Integer with the position from the start of the SNP of interest.
#' @export
drawAlignedMotifSNP <- function(mot.al.pos, length, posSNP) {
  motifStack::plotMotifLogoStack(mot.al.pos, xaxis=TRUE, yaxis=FALSE, ncex=0.9, ic.scale=TRUE, ylab=NA)
  grid.rect(x=1/length*posSNP, y=0.994, width=1/length, height=0.99, just=c("left", "top"),
            gp=gpar(lwd=2, lty="dashed", col=NA, fill="dark red", alpha=0.1))
}
