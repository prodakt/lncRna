#' BestTool.comb function
#'
#' This function computes confusion matrices and accuracy values for all combinations of all selected tools.
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @keywords confiusion matrix accuracy precision lncRNA
#' @export
#' @examples
#' BP.cmb <- BestTool.comb(BestPat = BestPat3)
#' BP.cmb

BestTool.comb <- function(BestPat, selectComb){

  '%!in%' <- function(x,y)!('%in%'(x,y))

  dat <- as.factor(BestPat[,selectComb[1]])
  levels(dat) <- c("0","1")
  ref <- as.factor(BestPat$isNC)
  levels(ref) <- c("0","1")

  cm <- caret::confusionMatrix(data = dat,
                        reference = ref, mode = "prec_recall", positive = "1")

  BP.cmb <- data.frame(torem = cm$overall)

  for(i in which(colnames(BestPat) %in% selectComb)) {

    dat <- as.factor(BestPat[,i])
    ref <- as.factor(BestPat$isNC)
    levels(dat) <- c("0","1")
    levels(ref) <- c("0","1")

    cm <- caret::confusionMatrix(data = dat,
                          reference = ref,
                          mode = "prec_recall", positive = "1")
    cm <- data.frame(cm$overall)
    cm <- round(cm, 4)

    BP.cmb[,ncol(BP.cmb) +1] <- cm
    names(BP.cmb)[ncol(BP.cmb)] <- colnames(BestPat)[i]
  }
  BP.cmb <- BP.cmb[,colnames(BP.cmb) %!in% "torem"]
  return(BP.cmb)
}
