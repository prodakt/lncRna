# --- Helper Functions ---

#' Interactively Select Tools or Default to All in Non-interactive Mode
#'
#' Prompts the user to select tools from a list if the session is interactive.
#' In a non-interactive session, it defaults to returning all available tools.
#'
#' @param availableTools A character vector of available tool names.
#' @return A character vector of selected tool names. Returns all tools in
#'   non-interactive mode. Returns `NULL` if interactive selection is cancelled.
#' @keywords internal
#' @noRd
selectToolsInteractively <- function(availableTools) {
  if (!interactive()) {
    # Non-interactive mode (e.g., R CMD check, vignettes): return all tools.
    return(availableTools)
  }

  # --- Interactive part ---
  message("Available tools:")
  toolOptions <- c("All Tools", availableTools)
  choicesText <- paste0(seq_along(toolOptions), ": ", toolOptions)
  message(paste(choicesText, collapse = "\n"))

  prompt <- "Enter numbers for tools to analyze (e.g., 1 for All, or 2,4): "
  userInput <- readline(prompt = prompt)
  userInput <- gsub("\\s+", "", userInput)

  if (userInput == "") {
    warning("No selection made.", call. = FALSE)
    return(NULL)
  }

  selectedIndices <- suppressWarnings(as.numeric(strsplit(userInput, ",")[[1]]))
  selectedIndices <- selectedIndices[!is.na(selectedIndices)]

  if (length(selectedIndices) == 0) {
    warning("Invalid selection. Please enter numbers.", call. = FALSE)
    return(NULL)
  }

  if (1 %in% selectedIndices) {
    return(availableTools)
  }

  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(toolOptions)]
  if (length(validIndices) == 0) {
    warning("Selected numbers do not correspond to available tools.", call. = FALSE)
    return(NULL)
  }
  return(availableTools[validIndices - 1])
}

#' Validate Provided Tool Names
#'
#' Checks if provided tool names exist in a list of available tools.
#'
#' @param toolsToValidate Character vector of tool names to check.
#' @param availableTools Character vector of all available tool names.
#' @return A character vector of validated tool names. Issues a warning for
#'   any names that were not found.
#' @keywords internal
#' @noRd
validateToolNames <- function(toolsToValidate, availableTools) {
  invalidTools <- setdiff(toolsToValidate, availableTools)
  if (length(invalidTools) > 0) {
    warning(
      "The following tools were not found and will be skipped: ",
      paste(invalidTools, collapse = ", "), call. = FALSE
    )
  }
  validatedTools <- intersect(toolsToValidate, availableTools)
  return(validatedTools)
}

#' Validate Input List Structure and Basic Data
#'
#' Checks if an input is a list and contains required elements.
#'
#' @param inputList The list to validate.
#' @param requiredElements A character vector of names that must be present.
#' @return `TRUE` if validation passes, `FALSE` otherwise (with warnings).
#' @keywords internal
#' @noRd
validateInputList <- function(inputList, requiredElements) {
  if (!is.list(inputList)) {
    stop("Input must be a list.")
  }

  missingElements <- setdiff(requiredElements, names(inputList))
  if (length(missingElements) > 0) {
    warning(
      "Input list is missing required elements: ",
      paste(missingElements, collapse = ", "), call. = FALSE
    )
    return(FALSE)
  }

  # Check for non-empty core data element ('isNC' or 'seqIDs')
  idSource <- intersect(c("isNC", "seqIDs"), names(inputList))[1]
  if (is.na(idSource)) {
    # This case is unlikely if requiredElements is checked, but good practice
    warning("Cannot determine number of sequences (missing 'isNC' or 'seqIDs').", call. = FALSE)
    return(FALSE)
  }
  nSeqs <- length(inputList[[idSource]])
  if (nSeqs == 0) {
    warning("Input list contains no sequences based on '", idSource, "'.", call. = FALSE)
    return(FALSE)
  }


  # Check for presence and basic structure of tool predictions if required
  if ("tools" %in% requiredElements) {
    if (!is.list(inputList$tools)) {
      warning("The 'tools' element in the input list is not a list.", call. = FALSE)
      return(FALSE)
    }
    if (length(inputList$tools) > 0) {
      firstToolLen <- length(inputList$tools[[1]])
      if (firstToolLen != nSeqs) {
        warning(
          "Length of first tool vector (", firstToolLen,
          ") does not match number of sequences (", nSeqs,
          ") based on '", idSource, "'.", call. = FALSE
        )
        return(FALSE) # Mismatch in length is a critical failure
      }
    }
  }
  return(TRUE)
}



#' Calculate Statistics for a Single Tool
#'
#' Computes confusion matrix statistics for one tool using caret. Returns a named
#' vector containing all available metrics, including NAs.
#'
#' @param toolName The name of the tool to analyze.
#' @param predictions A numeric vector (0/1) of predictions for the tool.
#' @param referenceFactor A factor vector of true labels (levels "0", "1").
#' @return A named numeric vector containing combined overall and by-class statistics
#'         for the tool (may include NAs), or NULL if calculation fails catastrophically.
#' @keywords internal statistics helper caret confusion matrix
#' @importFrom caret confusionMatrix
#' @noRd
#'
calculateSingleToolStats <- function(toolName, predictions, referenceFactor) {
  # Check prediction validity
  if (is.null(predictions) || length(predictions) != length(referenceFactor)) {
    warning(paste("Skipping tool '", toolName, "' due to missing or mismatched predictions."), call. = FALSE)
    return(NULL)
  }
  if (!all(predictions %in% c(0, 1))) {
    warning(paste("Skipping tool '", toolName, "' as predictions are not all 0 or 1."), call. = FALSE)
    return(NULL)
  }
  
  dat_factor <- factor(predictions, levels = levels(referenceFactor)) # Use same levels as reference
  
  # Use tryCatch for robustness
  cm_result <- tryCatch({
    caret::confusionMatrix(data = dat_factor,
                           reference = referenceFactor,
                           mode = "everything", # Get all available stats
                           positive = "1") # Assume '1' (non-coding) is positive class
  }, error = function(e) {
    warning(paste("Could not compute confusion matrix for tool '", toolName, "'. Error: ", e$message), call. = FALSE)
    return(NULL) # Return NULL if confusionMatrix itself fails
  })
  
  if (is.null(cm_result)) {
    return(NULL)
  }
  
  # Extract and combine stats
  stats_overall <- cm_result$overall
  stats_byClass <- cm_result$byClass
  # Combine into a single named vector
  combined_stats_vector <- c(stats_overall, stats_byClass)
  
  # Ensure the output is numeric, converting potential non-numeric outputs to NA
  # This handles cases where caret might return characters like "NaN"
  numeric_stats_vector <- suppressWarnings(as.numeric(combined_stats_vector))
  names(numeric_stats_vector) <- names(combined_stats_vector)
  
  return(numeric_stats_vector) # Return the named vector (may contain NAs)
}