#' SumSingleTools function
#'
#' This function summarise coding potential results for single tools.
#' @param CodPot.tbl2 the table with categorical values (0,1) where 0 means high coding potential and 1 means non-coding sequence.
#' @param nc_test the list of names of non-coding sequences for accuracy test (should be as the part of sequence names in the CodPot.tbl2)
#' @param cds_test the list of names of protein coding sequences for accuracy test (should be as the part of sequence names in the CodPot.tbl2)
#' @keywords training test lncRNA
#' @export
#' @examples
#' BestPat1 <- SumSingleTools(CodPot.tbl2 = tbl2, nc_test = nc_tt$nc.test, cds_test = cds_tt$cds.test)

SumSingleTools <- function(CodPot.tbl2, nc_test, cds_test){
  BestPat <- CodPot.tbl2
  usedTools <- colnames(CodPot.tbl2)[-1]
  ntools <- ncol(BestPat)-1
  BestPat$type <- "other"
  BestPat[BestPat$seqIDs %in% nc_test,]$type <- "nc"
  BestPat[BestPat$seqIDs %in% cds_test,]$type <- "cds"
  BestPat$isNC <- NA
  BestPat[BestPat$type %in% "cds",]$isNC <- 0
  BestPat[BestPat$type %in% "nc",]$isNC <- 1
  BestPat$sums <- rowSums(BestPat[,2:(ntools+1)])
  # removes those that are neither nc nor cds
  BestPat <- BestPat[!is.na(BestPat$isNC),]
  return(BestPat)
}
