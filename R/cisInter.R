#' cisInter function
#' This function reads several outputs from coding potential programs and join all tesults in table, where 0 is coding and 1 - noncoding
#' @param lncRNAs the list of lncRNA IDs (genes or transcripts)
#' @param mRNAs the list of mRNA IDs (genes or transcripts)
#' @param FEELnc.classes the "classes" output file from FEELnc
#' @param is.best filtering by the first collumn "is.best" (see https://github.com/tderrien/FEELnc)
#' @param lncRNA.level choose the "gene" or "transcript" level for lncRNA filtering
#' @param mRNA.level choose the "gene" or "transcript" level for mRNA filtering
#' @param max.dist the distance for cis-relation (default value is 100 000)
#' @keywords cis interactions lncRNA
#' @export
#' @examples
#' cisInter()
#'
cisInter <-function(lncRNAs = NULL, mRNAs = NULL, FEELnc.classes, is.best=T, lncRNA.level = "transcript", mRNA.level = "gene", max.dist = 100000){
  cis <- read.table(FEELnc.classes, header = T)
  if (is.best) cis <- cis[cis$isBest %in% 1,]

  if (length(lncRNAs) > 0) {
    if (lncRNA.level == "transcript") {
      cis <- cis[cis$lncRNA_transcript %in% lncRNAs,]
    } else {
      cis <- cis[cis$lncRNA_gene %in% lncRNAs,]
    }

  }

  if (length(mRNAs) > 0) {
    if (mRNA.level == "transcript") {
      cis <- cis[cis$partnerRNA_transcript %in% mRNAs,]
    } else {
      cis <- cis[cis$partnerRNA_gene %in% mRNAs,]
    }

  }
  cis <- cis[cis$distance <= max.dist,]
  names(cis)[names(cis) == 'lncRNA_gene'] <- 'lncRNA.id'
  names(cis)[names(cis) == 'partnerRNA_gene'] <- 'targetRNA.id'
  return(cis)
}
