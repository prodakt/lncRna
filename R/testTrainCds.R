#' Split Protein-Coding Sequences into Test and Train Sets
#'
#' Splits a set of protein-coding sequence identifiers into training and
#' testing subsets based on a specified proportion. For reproducible splits,
#' set the random seed using `set.seed()` before calling this function.
#'
#' @param cdsIds A character vector of protein-coding sequence identifiers, or
#'   an object (like a list from `seqinr::read.fasta`) from which names can be
#'   extracted.
#' @param trainProportion The proportion (0 to 1) of sequences to allocate to the
#'   training set (default: 0.6).
#'
#' @return A list containing two character vectors:
#'   \item{train}{Sequence identifiers allocated to the training set.}
#'   \item{test}{Sequence identifiers allocated to the testing set.}
#' @keywords split data partitioning sampling machine learning
#' @export
#' @examples
#' # --- 1. Define a character vector of sequence IDs ---
#' cdsIdsVec <- paste0("cds_seq_", 1:200)
#'
#' # --- 2. Run the function with a 50/50 split ---
#' set.seed(456) # for a reproducible example
#' cdsSplit <- splitCdsTestTrain(cdsIds = cdsIdsVec, trainProportion = 0.5)
#'
#' # --- 3. Verify the output ---
#' print(paste("Training set size:", length(cdsSplit$train)))
#' print(paste("Test set size:", length(cdsSplit$test)))
#' head(cdsSplit$train)
#'
splitCdsTestTrain <- function(cdsIds, trainProportion = 0.6) {
  # --- 1. Input Validation and Name Extraction ---
  if (is.character(cdsIds)) {
    sequenceNames <- cdsIds
  } else if (!is.null(names(cdsIds))) {
    sequenceNames <- names(cdsIds)
  } else {
    stop("'cdsIds' must be a character vector or an object with names.")
  }

  if (length(sequenceNames) == 0) {
    warning("No sequence identifiers found in 'cdsIds' input.")
    return(list(train = character(0), test = character(0)))
  }

  stopifnot(
    "'trainProportion' must be a number between 0 and 1" =
      (is.numeric(trainProportion) && trainProportion >= 0 && trainProportion <= 1)
  )

  # --- 2. Splitting Logic ---
  nTotal <- length(sequenceNames)
  nTrain <- round(nTotal * trainProportion)

  if (nTrain <= 0) {
    trainIds <- character(0)
    testIds <- sequenceNames
  } else if (nTrain >= nTotal) {
    trainIds <- sequenceNames
    testIds <- character(0)
  } else {
    trainIndices <- sample.int(nTotal, nTrain)
    trainIds <- sequenceNames[trainIndices]
    testIds <- sequenceNames[-trainIndices]
  }

  # --- 3. Return the list with standardized names ---
  return(list(train = trainIds, test = testIds))
}
