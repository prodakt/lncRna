#'
#' This function allows you to extract some statistics from reference GTF file
#' @param refGTF is the reference GTF imported with 'importGFF()' function from 'rtracklayer' package
#' @param genes can have to values TRUE - extract information about genes and FALSE - extract informations about transcripts
#' @keywords GTF lncRNA
#' @export
#' @examples
#' refBiotypes('GTF from stringtie by 'importGFF()')


refBiotypes <- function(refGTF, genes = F){
  known_biotypes <- as.data.frame(refGTF)

  if (genes) {
    known_biotypes <- unique(known_biotypes[,c("gene_id", "gene_biotype")])
    known_biotypes <- known_biotypes[!is.na(known_biotypes$gene_id),]
  } else {
    known_biotypes <- unique(known_biotypes[,c("transcript_id", "transcript_biotype")])
    known_biotypes <- known_biotypes[!is.na(known_biotypes$transcript_id),]
  }

  return(known_biotypes)
}

