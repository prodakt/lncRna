#' BestTool.comb function
#'
#' This function computes confusion matrices and accuracy values for all combinations of all selected tools.
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @keywords confiusion matrix accuracy precision lncRNA
#' @export
#' @examples
#' BP.cmb <- BestTool.comb(BestPat = BestPat3)
#' BP.cmb
#'

BestTool.comb <- function(BestPat){
  cm <- confusionMatrix(data = as.factor(BestPat[,colnames(cpAllComb)[1]]),
                        reference = as.factor(BestPat$isNC), mode = "prec_recall")
  BP.cmb <- data.frame(torem = cm$overall)

  for(i in which(colnames(BestPat) %in% colnames(cpAllComb))) {
    cm <- confusionMatrix(data = as.factor(BestPat[,i]),
                          reference = as.factor(BestPat$isNC),
                          mode = "prec_recall")
    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)

    BP.cmb[,ncol(BP.cmb) +1] <- cm
    names(BP.cmb)[ncol(BP.cmb)] <- colnames(BestPat)[i]
  }
  BP.cmb <- BP.cmb[,colnames(BP.cmb) %!in% "torem"]
  return(BP.cmb)
}
