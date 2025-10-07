#' Split Non-Coding Sequences into Test and Train Sets
#'
#' Splits a set of non-coding sequence identifiers into training and testing
#' subsets based on a specified proportion. For reproducible splits, set the
#' random seed using `set.seed()` before calling this function.
#'
#' @param ncIds A character vector of non-coding sequence identifiers, or an
#'   object (like a list from `seqinr::read.fasta`) from which names can be
#'   extracted.
#' @param trainProportion The proportion (0 to 1) of sequences to allocate to the
#'   training set (default: 0.6).
#'
#' @return A list containing two character vectors:
#'   \item{train}{Sequence identifiers allocated to the training set.}
#'   \item{test}{Sequence identifiers allocated to the testing set.}
#' @keywords split data partitioning sampling machine learning lncRNA
#' @export
#' @examples
#' # --- 1. Define a character vector of sequence IDs ---
#' ncIdsVec <- paste0("nc_seq_", 1:100)
#'
#' # --- 2. Run the function with a 70/30 split ---
#' set.seed(123) # for a reproducible example
#' ncSplit <- splitNcTestTrain(ncIds = ncIdsVec, trainProportion = 0.7)
#'
#' # --- 3. Verify the output ---
#' print(paste("Training set size:", length(ncSplit$train)))
#' print(paste("Test set size:", length(ncSplit$test)))
#' head(ncSplit$train)
#'
splitNcTestTrain <- function(ncIds, trainProportion = 0.6) {
  # --- 1. Input Validation and Name Extraction ---
  if (is.character(ncIds)) {
    sequenceNames <- ncIds
  } else if (!is.null(names(ncIds))) {
    sequenceNames <- names(ncIds)
  } else {
    stop("'ncIds' must be a character vector or an object with names.")
  }

  if (length(sequenceNames) == 0) {
    warning("No sequence identifiers found in 'ncIds' input.")
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
