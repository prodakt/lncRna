#' BestTool.atleast function
#'
#' This function computes confusion matrices and accuracy values for all variants of "at least n tools".
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @keywords confiusion matrix accuracy precision lncRNA
#' @export
#' @examples
#' BP.atl <- BestTool.atleast(BestPat = BestPat2)
#' BP.atl

BestTool.atleast <- function(BestPat){

  dat <- as.factor(BestPat[,2])
  levels(dat) <- c("0","1")
  ref <- as.factor(BestPat$isNC)
  levels(ref) <- c("0","1")

  cm <- caret::confusionMatrix(data = dat,
                        reference = ref, mode = "prec_recall", positive = "1")

  BP.atl <- data.frame(torem = cm$overall)

  for(i in which(startsWith(colnames(BestPat), "atl"))) {
    dat <- as.factor(BestPat[,i])
    ref <- as.factor(BestPat$isNC)
    levels(dat) <- c("0","1")
    levels(ref) <- c("0","1")

    cm <- caret::confusionMatrix(data = dat,
                          reference = ref,
                          mode = "prec_recall", positive = "1")

    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)
    BP.atl[,ncol(BP.atl) +1] <- cm
    names(BP.atl)[ncol(BP.atl)] <- colnames(BestPat)[i]
  }
  BP.atl <- BP.atl[,colnames(BP.atl) %!in% "torem"]

  return(BP.atl)
}
