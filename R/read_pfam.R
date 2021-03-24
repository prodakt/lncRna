#' A read.pfam function
#'
#' This function reads hmmscan pfam output and list all noncoding transcripts IDs
#' @param pfam_outfile is the hmmscan pfam output file localization including filename
#' @param eval_cutoff is the e-value cutoff
#' @keywords pfam lncRNA
#' @export
#' @examples
#' read.pfam()
#'
read.pfam <- function(pfam_outfile, eval_cutoff = 10e-3){
  require('rhmmer')
  pfam <- read_domtblout(file = pfam_outfile)
  pfam <- pfam[pfam$sequence_evalue < eval_cutoff,]$query_name
  pfam <- gsub("\\.p.*$","",pfam)
  pfam <- unique(pfam)
return(pfam)
}


