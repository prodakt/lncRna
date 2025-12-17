#' Estimate Trans-Acting Interactions Based on Expression Correlation
#'
#' This function estimates trans-interactions between lncRNAs and potential
#' target genes by calculating expression correlation. It uses the `Hmisc`
#' package for efficient correlation and p-value computation.
#'
#' @param exprMatrix A numeric matrix of expression values, where rows are genes
#'   (or transcripts) and columns are samples. Rownames must be gene identifiers.
#' @param corMethod The correlation method to use, either "pearson" (default) or
#'   "spearman".
#' @param rValueCutoff The absolute correlation coefficient cutoff. Interactions
#'   with an |r-value| greater than or equal to this are retained. Defaults to 0.7.
#' @param pValueCutoff The p-value cutoff. Interactions with a p-value less than
#'   this are retained. Defaults to 0.05.
#' @param lncRnaIds A character vector of lncRNA identifiers. These must match
#'   rownames in `exprMatrix`. If `NULL` (default), all genes are considered
#'   potential lncRNAs.
#' @param targetRnaIds A character vector of target gene identifiers. These must
#'   match rownames in `exprMatrix`. If `NULL` (default), all genes are
#'   considered potential targets.
#' @param fullCorrMatrixFile An optional character string specifying a file path
#'   to save the complete, unfiltered correlation matrix.
#'
#' @return A data.frame containing significant trans-interactions with columns:
#'   `lncRNAId`, `targetRNAId`, `rValue`, and `pValue`.
#'
#' @importFrom Hmisc rcorr
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @export
#' @examples
#' # --- 1. Create a sample expression matrix for a reproducible example ---
#' set.seed(123)
#' sampleExprMatrix <- matrix(
#'   rnorm(50), nrow = 5,
#'   dimnames = list(
#'     c("LNC_1", "LNC_2", "GENE_A", "GENE_B", "GENE_C"), # Rownames are gene IDs
#'     paste0("Sample_", 1:10)                             # Colnames are sample IDs
#'   )
#' )
#' # Introduce a strong correlation for demonstration
#' sampleExprMatrix["LNC_1",] <- sampleExprMatrix["GENE_A",] + rnorm(10, 0, 0.1)
#'
#' # --- 2. Define lncRNA and target gene lists ---
#' lncIds <- c("LNC_1", "LNC_2")
#' targetIds <- c("GENE_A", "GENE_B", "GENE_C")
#'
#' # --- 3. Run the analysis ---
#' if (requireNamespace("Hmisc", quietly = TRUE) &&
#'     requireNamespace("reshape2", quietly = TRUE)) {
#'   transInteractions <- calculateTransInteractions(
#'     exprMatrix = sampleExprMatrix,
#'     rValueCutoff = 0.9,
#'     lncRnaIds = lncIds,
#'     targetRnaIds = targetIds
#'   )
#'   print(transInteractions)
#'
#'   # --- 4. Example of saving the full correlation matrix to a temporary file ---
#'   tempFile <- tempfile(fileext = ".tsv")
#'   calculateTransInteractions(
#'     exprMatrix = sampleExprMatrix,
#'     lncRnaIds = lncIds,
#'     targetRnaIds = targetIds,
#'     fullCorrMatrixFile = tempFile
#'   )
#'   # In a real session, you could inspect this file.
#'   unlink(tempFile) # Clean up
#' }
calculateTransInteractions <- function(exprMatrix, corMethod = "pearson",
                                       rValueCutoff = 0.7, pValueCutoff = 0.05,
                                       lncRnaIds = NULL, targetRnaIds = NULL,
                                       fullCorrMatrixFile = NULL) {

  # --- 1. Input Validation ---
  corMethod <- match.arg(corMethod, c("pearson", "spearman"))
  stopifnot(
    "'exprMatrix' must be a numeric matrix with rownames" =
      (is.matrix(exprMatrix) && is.numeric(exprMatrix) && !is.null(rownames(exprMatrix)))
  )

  allIdsInMatrix <- rownames(exprMatrix)
  if (is.null(lncRnaIds)) {
    lncRnaIds <- allIdsInMatrix
  }
  if (is.null(targetRnaIds)) {
    targetRnaIds <- allIdsInMatrix
  }

  stopifnot(
    "All 'lncRnaIds' must be present in 'exprMatrix' rownames" = all(lncRnaIds %in% allIdsInMatrix),
    "All 'targetRnaIds' must be present in 'exprMatrix' rownames" = all(targetRnaIds %in% allIdsInMatrix)
  )

  # --- 2. Prepare Data and Calculate Correlations ---
  idsToKeep <- unique(c(as.character(lncRnaIds), as.character(targetRnaIds)))
  subExprMatrix <- exprMatrix[idsToKeep, , drop = FALSE]

  # Hmisc::rcorr requires variables in columns, so we transpose the matrix
  corrResult <- Hmisc::rcorr(t(subExprMatrix), type = corMethod)

  if (!is.null(fullCorrMatrixFile)) {
    utils::write.table(as.data.frame(corrResult$r), fullCorrMatrixFile, quote = FALSE)
  }

  # --- 3. Process and Filter Results ---
  # Melt correlation and p-value matrices into long format
  transRes <- reshape2::melt(corrResult$r, varnames = c("Var1", "Var2"), value.name = "rValue")
  pValues <- reshape2::melt(corrResult$P, varnames = c("Var1", "Var2"), value.name = "pValue")
  transRes$pValue <- pValues$pValue

  # Filter for lncRNA-target pairs
  transRes <- transRes[transRes$Var1 %in% lncRnaIds & transRes$Var2 %in% targetRnaIds, ]

  # Filter by p-value and r-value cutoffs
  transRes <- transRes[transRes$pValue < pValueCutoff, ]
  transRes <- transRes[abs(transRes$rValue) >= rValueCutoff, ]

  # Remove self-correlations if an ID is in both lists
  transRes <- transRes[as.character(transRes$Var1) != as.character(transRes$Var2), ]

  # Standardize column names
  colnames(transRes)[1:2] <- c("lncRNAId", "targetRNAId")

  message(paste0("Number of trans-interactions found: ", nrow(transRes)))

  return(transRes)
}
