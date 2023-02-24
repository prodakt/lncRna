#' BestTool.atleast function
#'
#' This function computes confusion matrices and accuracy values for all variants of "at least n tools".
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @keywords confiusion matrix accuracy precision lncRNA
#' @export
#' @examples
#' BP.atl <- BestTool.atleast(BestPat = BestPat2)
#' BP.atl
#'

BestTool.atleast <- function(BestPat){
  cm <- confusionMatrix(data = as.factor(BestPat[,2]),
                        reference = as.factor(BestPat$isNC), mode = "prec_recall")
  BP.atl <- data.frame(torem = cm$overall)

  for(i in which(startsWith(colnames(BestPat), "atl"))) {
    cm <- confusionMatrix(data = as.factor(BestPat[,i]),
                          reference = as.factor(BestPat$isNC),
                          mode = "prec_recall")
    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)
    BP.atl[,ncol(BP.atl) +1] <- cm
    names(BP.atl)[ncol(BP.atl)] <- colnames(BestPat)[i]
  }
  BP.atl <- BP.atl[,colnames(BP.atl) %!in% "torem"]

  return(BP.atl)
}
