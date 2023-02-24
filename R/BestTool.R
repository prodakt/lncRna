#' BestTool function
#'
#' This function computes confusion matrix and accuracy for single methods.
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @param tools the list of names of tools to compute the accuracy
#' @keywords training test lncRNA
#' @export
#' @examples
#' selectedTools <- c("CPC2", "PLEK", "FEELnc", "CPAT", "CNCI", "LncFinder")
#' BP.cpt <- BestTool(BestPat = BestPat1, selectedTools)
#' BP.cpt

BestTool <- function(BestPat, tools){
  ntools <- length(tools)
  cm <- confusionMatrix(data = as.factor(BestPat[,2]),
                        reference = as.factor(BestPat$isNC), mode = "prec_recall")
  BP.cpt <- data.frame(torem = cm$overall)
  for(i in which(colnames(BestPat) %in% usedTools)) {
    cm <- confusionMatrix(data = as.factor(BestPat[,i]),
                          reference = as.factor(BestPat$isNC),
                          mode = "prec_recall")
    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)
    BP.cpt <- cbind(BP.cpt,cm)
    #  colnames(BP.cpt)[ncol(BP.cpt)] <- cat(colnames(dane)[i])
  }
  BP.cpt <- BP.cpt[,colnames(BP.cpt) %!in% "torem"]
  colnames(BP.cpt)[1:ntools] <- colnames(BestPat)[c(1:ntools+1)]
  return(BP.cpt)
}
