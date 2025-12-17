#' Find Potential trans-Regulatory Interactions
#'
#' Estimates trans-interactions between lncRNAs and target RNAs based on
#' expression correlation. This function requires the `Hmisc` and `reshape2`
#' packages.
#'
#' @param exprMatrix A numeric matrix or data.frame of expression values.
#'   Rownames should contain gene/transcript IDs and columns should be samples.
#' @param corMethod Correlation method: `"pearson"` (default) or `"spearman"`.
#' @param rval The cutoff for the correlation coefficient (default: 0.7).
#' @param pval The cutoff for the p-value (default: 0.05).
#' @param lncRnaList A list of lncRNA gene/transcript IDs. Must be present in
#'   `rownames(exprMatrix)`. If `NULL`, all rownames are considered.
#' @param tarRnaList A list of target gene/transcript IDs. Must be present in
#'   `rownames(exprMatrix)`. If `NULL`, all rownames are considered.
#' @param fullCorMatrixFile An optional file path to save the full correlation matrix.
#'
#' @return A `data.frame` of significant trans-interactions with columns for
#'   lncRNA ID, target RNA ID, r-value, and p-value.
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#'
#' @examples
#' # --- 1. Create a mock expression matrix ---
#' set.seed(123)
#' lnc_genes <- paste0("LNC", 1:5)
#' target_genes <- paste0("TARGET", 1:20)
#' all_genes <- c(lnc_genes, target_genes)
#' mockExprMatrix <- matrix(rnorm(25 * 10), nrow = 25, ncol = 10,
#'                          dimnames = list(all_genes, paste0("Sample", 1:10)))
#' mockExprMatrix["LNC1", ] <- mockExprMatrix["TARGET1", ] * 2 + rnorm(10, 0, 0.2)
#'
#' # --- 2. Run the function ---
#' trans_interactions <- findTransInteractions(
#'   exprMatrix = mockExprMatrix,
#'   lncRnaList = lnc_genes,
#'   tarRnaList = target_genes,
#'   rval = 0.9,
#'   pval = 0.05
#' )
#' print(trans_interactions)
#'
findTransInteractions <- function(exprMatrix, corMethod = "pearson", rval = 0.7,
                                  pval = 0.05, lncRnaList = NULL,
                                  tarRnaList = NULL, fullCorMatrixFile = NULL) {

    if (!is.matrix(exprMatrix) && !is.data.frame(exprMatrix)) {
        stop("'exprMatrix' must be a matrix or a data.frame.")
    }

    if (is.null(lncRnaList)) {
        lncRnaList <- rownames(exprMatrix)
    }
    if (is.null(tarRnaList)) {
        tarRnaList <- rownames(exprMatrix)
    }

    exprMatrixFiltered <- exprMatrix[rownames(exprMatrix) %in% c(as.character(lncRnaList), as.character(tarRnaList)), ]

    sim_matrix <- Hmisc::rcorr(t(as.matrix(exprMatrixFiltered)), type = corMethod)

    if (!is.null(fullCorMatrixFile)) {
        utils::write.table(as.data.frame(sim_matrix$r), fullCorMatrixFile, quote = FALSE)
    }

    trans_res <- reshape2::melt(sim_matrix$r)
    p <- reshape2::melt(sim_matrix$P)

    trans_res$pValue <- p$value

    colnames(trans_res) <- c("lncRNAId", "targetRNAId", "rValue", "pValue")

    trans_res <- trans_res[trans_res$lncRNAId %in% lncRnaList, ]
    trans_res <- trans_res[trans_res$targetRNAId %in% tarRnaList, ]

    trans_res <- trans_res[which(trans_res$pValue < pval), ]
    trans_res <- trans_res[which(abs(trans_res$rValue) >= rval), ]

    trans_res <- trans_res[!is.na(trans_res$lncRNAId), ]

    message("Number of trans lncRNA vs. target genes interactions: ", nrow(trans_res))

    return(trans_res)
}
