
#' A TransAct function
#'
#' This function allows to estimate trans interactions between lncRNA's and potential target genes basing on the expression/coexpression levels.
#' This function requires the "Hmisc" library to be installed
#' @param expr.matrix espression matrix; rownames should contains the genes or transcripts names and collumns shoud be the samples.
#' @param cor.method correlation metthod: "pearson" or "spearman"
#' @param rval the cutoff of correlation coeficient; default value is 0.7
#' @param pval the cutoff of estimated p-value for any interaction; default value is 0.05
#' @param lncRNA.list the list of lncRNA genes/transcripts. The manes shoud be a part of "expr.matrix" rownames.
#' @param tarRNA.list the list of target genes/transcripts. The manes shoud be a part of "expr.matrix" rownames.
#' @param full.cor.matrix.filename
#' @keywords trans acting lncRNA
#' @export
#' @examples transNjJd <- TransAct(expr.matrix = em_NjJd, rval = 0.9, lncRNA.list = DELsNjJd, tarRNA.list = pcDEGsNjJd)
#' @examples head(transNjJd)
#' TransAct()
#'

TransAct <- function(expr.matrix, cor.method = "pearson", rval = 0.7, pval = 0.05, lncRNA.list=NULL, tarRNA.list=NULL, full.cor.matrix.filename = NULL){

  require("Hmisc")

  if (is.null(lncRNA.list)) {
    lncRNA.list=rownames(expr.matrix)
  }
  if (is.null(tarRNA.list)) {
    tarRNA.list=rownames(expr.matrix)
  }
  expr.matrix <- expr.matrix[rownames(expr.matrix) %in% c(as.character(lncRNA.list), as.character(tarRNA.list)),]
  sim_matrix <- rcorr(t(as.matrix(expr.matrix)), type = cor.method)
  if (!is.null(full.cor.matrix.filename)) {
    write.table(as.data.frame(sim_matrix$r), full.cor.matrix.filename, quote = F)
  }
  trans_res <- melt(sim_matrix$r)
  p <- melt(sim_matrix$P)
  trans_res$p.value <- p$value
  colnames(trans_res) <- c("lncRNA.id", "targetRNA.id", "r.value", "p.value")
  trans_res <- trans_res[trans_res$lncRNA.id %in% lncRNA.list,]
  trans_res <- trans_res[trans_res$targetRNA.id %in% tarRNA.list,]
  trans_res <- trans_res[trans_res$p.value < pval,]
  trans_res <- trans_res[abs(trans_res$r.value) >= rval,]

  message(paste0("Number of trans lncRNA vs. target genes interactions: ", nrow(trans_res)))

  return(trans_res)

}
