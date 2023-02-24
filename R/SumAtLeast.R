#' SumAtLeast function
#'
#' This function summarise coding potential results for several tools counting if any sequence is predicted by at least n tools.
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @param tools the list of names of tools to compute the accuracy
#' @keywords training test lncRNA
#' @export
#' @examples
#' BestPat2 <- SumAtLeast(BestPat = BestPat1, tools = selectedTools)
#' head(BestPat2)
#'
SumAtLeast <- function(BestPat, tools){
  ntools <- length(tools)
  for (n in 1:ntools) {
    BestPat[,ncol(BestPat) +1] <- ifelse(BestPat$sums >= n,1,0)
    names(BestPat)[ncol(BestPat)] <- paste0('atl', n)
  }
  return(BestPat)
}
