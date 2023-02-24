#' SumCombTools function
#'
#' This function summarise coding potential results for all combinations of selected tools.
#' @param BestPat the table with summarized data generated with SumSingleTools() function
#' @param tools the list of names of tools to compute the accuracy
#' @keywords training test lncRNA
#' @export
#' @examples
#' BestPat3 <- SumAtLeast(BestPat = BestPat2, tools = selectedTools)
#' head(BestPat3)
#'

SumCombTools <- function(BestPat, selectedTools){
  met2comb <- BestPat[,selectedTools]
  rownames(met2comb) <- BestPat$seqIDs

  n <- ncol(met2comb)
  l <- rep(list(0:1), n)
  combinations <- expand.grid(l)
  colnames(combinations) <- selectedTools

  # rownames of all combinations
  for (i in 1:nrow(combinations)) {
    rownames(combinations)[i] <- ifelse(sum(combinations[i,]) == 0, "zero",
                                        paste0(colnames(combinations)[as.logical(combinations[i,])], collapse = "+")
    )
  }
  combinations <- combinations[rowSums(combinations) > 1,] # remove combinations not  existing - with sums of "0" value

  # generate the table of predicted in combinations - intersections of venn
  cpAllComb <- data.frame(matrix(NA,    # Create empty data frame
                                 nrow = nrow(met2comb),
                                 ncol = nrow(combinations)))
  rownames(cpAllComb) <- rownames(met2comb)
  colnames(cpAllComb) <- rownames(combinations)

  for (col in 1:ncol(cpAllComb)) {
    cpAllComb[,col] = 0 # substiture "NAs" by "0" to know it was done
    name <- colnames(cpAllComb)[col] # take the column name
    tmp.tab <- met2comb[,as.logical(combinations[rownames(combinations) %in% name,])] # generate temporary table
    selected <- rownames(tmp.tab)[rowSums(tmp.tab) == ncol(tmp.tab)] # take the names of sequences, that meets the criterion of intersection
    cpAllComb[rownames(cpAllComb) %in% selected, name] <- 1 # ... and put "1" in these sequences
  }

  # compute cm for all combinations
  comb.names <- colnames(cpAllComb)
  BestPat <- merge(BestPat, cpAllComb, by.x = "seqIDs", by.y = 0, all.x = T)

  return(BestPat)
}
