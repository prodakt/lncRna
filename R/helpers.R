# --- Helper Functions (potentially shared with BestTool) ---

#' Interactively Select Tools
#'
#' Prompts the user to select tools from a list if the session is interactive.
#' "All Tools" is presented as the first option.
#'
#' @param available_tools A character vector of available tool names.
#' @return A character vector of selected tool names, or NULL if selection fails or is cancelled.
#' @keywords internal interactive helper selection
#'
.select_tools_interactively <- function(available_tools) {
  if (!interactive()) {
    # In non-interactive mode for checks, we can't prompt.
    # Return the first tool if available for example execution.
    warning("Non-interactive session detected in example. Selecting first available tool.", call. = FALSE)
    if(length(available_tools) > 0) {
      return(available_tools[1])
    } else {
      return(NULL)
    }
  }

  # --- Interactive part (will run only in interactive sessions) ---
  cat("Available tools:\n")
  # Put "All Tools" as option 1
  tool_options <- c("All Tools", available_tools)
  choice_indices <- seq_along(tool_options)
  print(paste(choice_indices, ": ", tool_options, sep = ""), quote = FALSE)

  prompt <- "Enter the numbers for the tools to analyze (e.g., 1 for All, or 2,4), then press Enter: "
  user_input <- readline(prompt = prompt)
  user_input <- gsub("\\s+", "", user_input) # Remove whitespace

  if (user_input == "") {
    warning("No selection made.", call. = FALSE)
    return(NULL)
  }

  selected_indices_str <- strsplit(user_input, ",")[[1]]
  selected_indices <- suppressWarnings(as.numeric(selected_indices_str))
  selected_indices <- selected_indices[!is.na(selected_indices)] # Remove NAs

  if (length(selected_indices) == 0) {
    warning("Invalid selection. Please enter numbers corresponding to the list.", call. = FALSE)
    return(NULL)
  }

  # Check if "All Tools" (index 1) was selected
  if (1 %in% selected_indices) {
    return(available_tools) # User selected "All Tools"
  } else {
    # User selected specific tools (indices 2 onwards map to available_tools)
    # Validate selected indices are within the range of actual tools (2 to length(tool_options))
    valid_indices <- selected_indices[selected_indices > 1 & selected_indices <= length(tool_options)]
    if (length(valid_indices) == 0) {
      warning("Selected numbers do not correspond to available tools.", call. = FALSE)
      return(NULL)
    }
    # Map indices back to available_tools (subtract 1)
    return(available_tools[valid_indices - 1])
  }
}

#' Validate Provided Tool Names
#'
#' Checks if the provided tool names exist in the list of available tools.
#'
#' @param tools_to_validate Character vector of tool names provided by the user.
#' @param available_tools Character vector of all available tool names.
#' @return Character vector of validated tool names that exist in available_tools.
#'         Issues warnings for invalid names.
#' @keywords internal validation helper
#'
.validate_tool_names <- function(tools_to_validate, available_tools) {
  invalid_tools <- setdiff(tools_to_validate, available_tools)
  if (length(invalid_tools) > 0) {
    warning(paste("The following provided tools were not found and will be skipped:",
                  paste(invalid_tools, collapse = ", ")), call. = FALSE)
  }
  validated_tools <- intersect(tools_to_validate, available_tools)
  return(validated_tools)
}

#' Validate Input List Structure and Basic Data
#'
#' Checks if the input is a list, contains required elements, and has sequence data.
#'
#' @param input_list The list to validate.
#' @param required_elements A character vector of names that must be present in the list.
#' @return TRUE if validation passes, FALSE otherwise (with warnings).
#' @keywords internal validation helper list structure
#'
.validate_input_list <- function(input_list, required_elements) {
  # Check if it's a list
  if (!is.list(input_list)) {
    stop("Input must be a list.") # Stop for fundamental type mismatch
  }
  # Check for required elements
  missing_elements <- setdiff(required_elements, names(input_list))
  if (length(missing_elements) > 0) {
    warning(paste("Input list is missing required elements:", paste(missing_elements, collapse = ", ")), call. = FALSE)
    return(FALSE)
  }

  # Determine sequence count (use seqIDs if present, otherwise isNC)
  n_seqs <- 0
  id_source <- NULL
  if ("seqIDs" %in% names(input_list)) {
    n_seqs <- length(input_list$seqIDs)
    id_source <- "seqIDs"
  } else if ("isNC" %in% names(input_list)) {
    n_seqs <- length(input_list$isNC)
    id_source <- "isNC"
  } else {
    warning("Cannot determine number of sequences (missing 'seqIDs' or 'isNC').", call. = FALSE)
    return(FALSE)
  }

  # Check if sequences are present
  if (n_seqs == 0) {
    warning(paste("Input list contains no sequences based on element:", id_source), call. = FALSE)
    return(FALSE)
  }

  # Check for presence and basic structure of tool predictions if 'tools' is required
  if ("tools" %in% required_elements) {
    if (!is.list(input_list$tools)) {
      warning("Input list$tools is not a list.", call. = FALSE)
      return(FALSE)
    }
    if (length(input_list$tools) == 0) {
      warning("Input list$tools contains no tool predictions.", call. = FALSE)
      # Depending on function, this might be acceptable, but usually not.
      # Consider if this should be FALSE or just a warning. For now, let it pass.
      # return(FALSE)
    } else {
      # Optional: Check if first tool vector matches sequence length
      # This helps catch major inconsistencies early.
      first_tool_len <- length(input_list$tools[[1]])
      if (first_tool_len != n_seqs) {
        warning(paste0("Length of first tool vector (", first_tool_len,
                       ") does not match number of sequences (", n_seqs,
                       ") based on ", id_source, "."), call. = FALSE)
        # Decide if this is a critical failure
        # return(FALSE)
      }
    }
  }

  return(TRUE) # Basic validation passed
}
