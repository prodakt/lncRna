#' Split Identifiers into Test and Train Sets
#'
#' Splits a set of identifiers (e.g., for genes, transcripts, or samples) into
#' training and testing subsets based on a specified proportion. For
#' reproducible splits, set the random seed using `set.seed()` before calling
#' this function.
#'
#' @param ids A character vector of identifiers, or an object (like a list from
#'   `seqinr::read.fasta`) from which names can be extracted.
#' @param trainProportion The proportion (0 to 1) of identifiers to allocate to
#'   the training set (default: 0.6).
#'
#' @return A list containing two character vectors:
#'   \item{train}{Identifiers allocated to the training set.}
#'   \item{test}{Identifiers allocated to the testing set.}
#' @keywords split data partitioning sampling machine learning
#' @export
#' @examples
#' # --- 1. Demonstrate splitting for non-coding sequences ---
#' ncIdsVec <- paste0("nc_seq_", 1:100)
#' set.seed(123) # for reproducibility
#' ncSplit <- splitTestTrain(ids = ncIdsVec, trainProportion = 0.7)
#' print("Non-coding split:")
#' print(paste("Training set size:", length(ncSplit$train)))
#' print(paste("Test set size:", length(ncSplit$test)))
#'
#' # --- 2. Demonstrate splitting for protein-coding sequences ---
#' cdsIdsVec <- paste0("cds_seq_", 1:200)
#' set.seed(456) # for reproducibility
#' cdsSplit <- splitTestTrain(ids = cdsIdsVec, trainProportion = 0.5)
#' print("Protein-coding split:")
#' print(paste("Training set size:", length(cdsSplit$train)))
#' print(paste("Test set size:", length(cdsSplit$test)))
#'
splitTestTrain <- function(ids, trainProportion = 0.6) {
  # --- 1. Input Validation and Name Extraction ---
  if (is.character(ids)) {
    sequenceNames <- ids
  } else if (!is.null(names(ids))) {
    sequenceNames <- names(ids)
  } else {
    stop("'ids' must be a character vector or an object with names.")
  }

  if (length(sequenceNames) == 0) {
    warning("No identifiers found in 'ids' input.")
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

  # --- 3. Return the list ---
  return(list(train = trainIds, test = testIds))
}
