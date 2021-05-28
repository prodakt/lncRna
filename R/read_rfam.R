#' A read.rfam function
#'
#' This function reads hmmscan rfam output and list all noncoding transcripts IDs
#' @param rfam_outfile is the hmmscan pfam output file localization including filename
#' @param eval_cutoff is the e-value cutoff
#' @keywords rfam lncRNA
#' @export
#' @examples
#' read.rfam()
#'
read.rfam <- function(rfam_outfile){
  require('rhmmer')
  rfam <- read_domtblout(file = rfam_outfile)
  #
  # rfam <- rfam[rfam$sequence_evalue < eval_cutoff,]$query_name
  # rfam <- rfam[rfam$sequence_evalue < eval_cutoff,]
  rfam <- as.data.frame(rfam)
  rfam <- unique(rfam)
  return(rfam)
}


