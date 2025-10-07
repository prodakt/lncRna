#' Summarize Predictions for Combinations of Tools
#'
#' Identifies sequences predicted as non-coding by specific combinations of
#' selected tools (i.e., the intersection of predictions). It considers
#' combinations of 2 or more tools.
#'
#' @param sumSingleToolsList A list from a precursor function, containing at
#'   least `$seqIDs` and `$tools` (a named list of 0/1 prediction vectors).
#' @param tools An optional character vector of tool names to use for
#'   combinations. If `NULL` (default) and the session is interactive, the user
#'   is prompted to select tools. In a non-interactive session, `NULL` results
#'   in all tools being used.
#'
#' @return A list containing:
#'   \item{seqIDs}{Original sequence identifiers.}
#'   \item{isNC}{Original numeric indicator (if present).}
#'   \item{type}{Original type annotation (if present).}
#'   \item{selectedToolsPredictions}{A sublist of prediction vectors for the
#'     selected tools.}
#'   \item{toolCombinations}{A list where each element is named after a tool
#'     combination (e.g., "ToolA+ToolB") and contains a binary vector (0/1)
#'     indicating a match for each sequence.}
#'   Returns `NULL` if validation fails or fewer than two tools are selected.
#' @keywords summary combination coding potential lncRNA intersection
#' @importFrom utils combn
#' @export
#' @examples
#' # --- 1. Create Example Data ---
#' set.seed(101)
#' nSeqExample <- 50
#' exampleInputList <- list(
#'   seqIDs = paste0("Seq", seq_len(nSeqExample)),
#'   tools = list(
#'     ToolX = sample(c(0, 1), nSeqExample, replace = TRUE),
#'     ToolY = sample(c(0, 1), nSeqExample, replace = TRUE),
#'     ToolZ = sample(c(0, 1), nSeqExample, replace = TRUE)
#'   ),
#'   type = sample(c("nc", "cds"), nSeqExample, replace = TRUE)
#' )
#' exampleInputList$isNC <- ifelse(exampleInputList$type == "nc", 1, 0)
#'
#' # --- 2. Run Non-interactive Example ---
#' # Analyze combinations of all available tools
#' resultsComb <- sumCombTools(sumSingleToolsList = exampleInputList)
#'
#' if (!is.null(resultsComb)) {
#'   print("Results for all tool combinations:")
#'   # Print the names of the generated combinations
#'   print(names(resultsComb$toolCombinations))
#'   # Print the head of a specific combination's result vector
#'   if ("ToolX+ToolY" %in% names(resultsComb$toolCombinations)) {
#'     print(head(resultsComb$toolCombinations[["ToolX+ToolY"]]))
#'   }
#' }
#'
#' # --- 3. Interactive Example ---
#' if (interactive()) {
#'   message("Running interactive example. You will be prompted for selection.")
#'   resultsInteractive <- sumCombTools(sumSingleToolsList = exampleInputList)
#' }
sumCombTools <- function(sumSingleToolsList, tools = NULL) {

  # --- 1. Input Validation ---
  requiredElements <- c("seqIDs", "tools")
  # Assuming validateInputList is a helper function in your package
  if (!validateInputList(sumSingleToolsList, requiredElements)) {
    return(NULL)
  }

  availableTools <- names(sumSingleToolsList$tools)
  if (is.null(availableTools) || length(availableTools) == 0) {
    warning("No tool predictions found in 'sumSingleToolsList$tools'.", call. = FALSE)
    return(NULL)
  }

  # --- 2. Tool Selection ---
  if (is.null(tools)) {
    # selectToolsInteractively should return all tools if non-interactive
    selectedToolNames <- selectToolsInteractively(availableTools)
  } else {
    selectedToolNames <- validateToolNames(tools, availableTools)
  }

  if (is.null(selectedToolNames) || length(selectedToolNames) < 2) {
    warning("At least two valid tools must be selected to form combinations.", call. = FALSE)
    return(NULL)
  }

  # --- 3. Prepare Data for Selected Tools ---
  selectedPredictionsList <- sumSingleToolsList$tools[selectedToolNames]
  nSeqs <- length(sumSingleToolsList$seqIDs)

  if (!all(sapply(selectedPredictionsList, length) == nSeqs)) {
    warning("One or more selected tools have a prediction vector of incorrect length.", call. = FALSE)
    return(NULL)
  }

  predictionsMatrix <- do.call(cbind, selectedPredictionsList)

  # --- 4. Generate Combinations and Calculate Matches ---
  toolCombinationsResults <- list()
  nSelected <- length(selectedToolNames)

  for (k in 2:nSelected) {
    toolIndicesCombinations <- utils::combn(seq_len(nSelected), k, simplify = FALSE)

    for (idxVector in toolIndicesCombinations) {
      currentCombinationNames <- selectedToolNames[idxVector]
      combinationId <- paste(currentCombinationNames, collapse = "+")

      subMatrix <- predictionsMatrix[, idxVector, drop = FALSE]
      isMatch <- rowSums(subMatrix) == k
      toolCombinationsResults[[combinationId]] <- as.integer(isMatch)
    }
  }

  # --- 5. Assemble Output List ---
  finalList <- list(
    seqIDs = sumSingleToolsList$seqIDs,
    isNC = sumSingleToolsList$isNC,
    type = sumSingleToolsList$type,
    selectedToolsPredictions = selectedPredictionsList,
    toolCombinations = toolCombinationsResults
  )

  return(finalList)
}
