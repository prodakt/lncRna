#' A filter.CisAct function
#'
#' This function reads several outputs from coding potential programs and join all tesults in table, where 0 is coding and 1 - noncoding
#' @param lncRNAs the list of lncRNA IDs (genes or transcripts)
#' @param mRNAs the list of mRNA IDs (genes or transcripts)
#' @param FEELnc.classes the "classes" output file from FEELnc
#' @param is.best filtering by the first collumn "is.best" (see https://github.com/tderrien/FEELnc#3--feelnc_classifierpl)
#' @param lncRNA.level choose the "gene" or "transcript" level for lncRNA filtering
#' @param mRNA.level choose the "gene" or "transcript" level for mRNA filtering
#' @param
#' @keywords cis acting lncRNA
#' @export
#' @examples
#' filter.CisAct()
#'
filter.CisAct <-function(lncRNAs = NULL, mRNAs = NULL, FEELnc.classes, is.best=T, lncRNA.level = "transcript", mRNA.level = "gene"){
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

  return(cis)
}
