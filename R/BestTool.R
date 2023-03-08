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
  dat <- as.factor(BestPat[,2])
  levels(dat) <- c("0","1")
  ref <- as.factor(BestPat$isNC)
  levels(ref) <- c("0","1")

  cm <- caret::confusionMatrix(data = dat,
                        reference = ref, mode = "prec_recall", positive = "1")

 BP.cpt <- data.frame(torem = cm$overall)

 #  BP.cpt <- rbind(data.frame(torem = cm$overall), data.frame(torem = cm$byClass))

  for(i in which(colnames(BestPat) %in% tools)) {
    dat <- as.factor(BestPat[,i])
    ref <- as.factor(BestPat$isNC)
    levels(dat) <- c("0","1")
    levels(ref) <- c("0","1")

    cm <- caret::confusionMatrix(data = dat,
                          reference = ref,
                          mode = "prec_recall", positive = "1")
    # cm <- rbind(data.frame(tmp = cm$overall), data.frame(tmp = cm$byClass))
    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)
    BP.cpt <- cbind(BP.cpt,cm)
  }
  BP.cpt <- BP.cpt[,colnames(BP.cpt) %!in% "torem"]
  colnames(BP.cpt)[1:ntools] <- colnames(BestPat)[c(1:ntools+1)]
  return(BP.cpt)
}
