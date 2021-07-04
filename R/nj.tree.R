#' Draw ggtree based phylogenetic tree using DNA or protein sequences
#'
#' @param file a file path. The file containg the sequences in fasta format.
#' @param align.method method for multiple sequences alignment, one of "muscle" and "DECIPHER"
#' @param html a boolean value to indicate if to show the alignment result in the explore
#' @param plot a boolean value if draw the tree using ggtree package
#'
#' @author Zhonghui Gai
#' @return a ggtree plot
#' @export
#'
#' @examples
#' nj.tree("prob.16S.fa")
nj.tree <- function(file,
                    align.method = "muscle",
                    html = TRUE,
                    plot = TRUE){
  seq <- Biostrings::readDNAStringSet(filepath = file, format = "fasta",
                                      skip = 0L, seek.first.rec = FALSE,
                                      use.names = TRUE)
  if (align.method == "muscle") {
    aln <- muscle::muscle(seq)
    dist <- Biostrings::stringDist(as(aln,"AAStringSet"), method="hamming")
  }else if (align.method == "DECIPHER") {
    aln <- DECIPHER::AlignSeqs(seq)
    Biostrings::writeXStringSet(aln, filepath = "temp.fa")
    dist <- seqinr::dist.alignment(x = seqinr::read.alignment(file = "temp.fa",
                                                              format = "fasta"),
                                   matrix = "similarity")
    if (TRUE) {file.remove("temp.fa")}
    if (html) {
      DECIPHER::BrowseSeqs(aln, highlight = 0)
    }
  }
  tree <- ape::nj(dist)
  if (plot) {
    p <- ggtree::ggtree(tree, layout = "ellipse", branch.length = "none") +
      ggtree::geom_tiplab() + ggtree::theme_tree2()
  }else{
    return(tree)
  }
  return(p)
}
