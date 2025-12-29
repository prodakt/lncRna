#' Handle Tool Selection and Validation (Internal Helper)
#'
#' A reusable helper function to manage the selection of tools for analysis.
#' It supports both user-provided lists and an interactive selection mode.
#'
#' @param availableTools A character vector of all tool names available for
#'   selection.
#' @param selectedTools A character vector of tool names provided by the user,
#'   or `NULL` to trigger default/interactive mode.
#'
#' @return A character vector of validated tool names to be used in an
#'   analysis. Returns an empty character vector (`character(0)`) if no
#'   valid tools are selected.
#' @noRd
handleToolSelection <- function(availableTools, selectedTools = NULL) {
    
    final_selected_tools <- character(0)
    
    if (is.null(selectedTools)) {
        if (interactive()) {
            message("Available tools:\n")
            tool_menu <- paste0(seq_along(availableTools), ": ", availableTools)
            message(paste(tool_menu, collapse = "\n"), "\n")
            
            input <- readline(
                "Enter tool numbers to analyze (e.g., '1 3'), or press Enter for all: "
            )
            
            if (input == "") {
                final_selected_tools <- availableTools
                message("All tools selected.\n")
            } else {
                # suppressWarnings() is used here to silence expected warnings 
                # (during as.integer(strsplit()) calls . These warnings are 
                # expected and do not impact the downstream analysis. This approach 
                # is chosen to maintain code readability and simplicity.
                indices <- suppressWarnings(as.integer(strsplit(input, "[ ,]+")[[1]]))
                valid_indices <- indices[!is.na(indices) & indices > 0 &
                                             indices <= length(availableTools)]
                final_selected_tools <- availableTools[valid_indices]
            }
        } else {
            # Default non-interactive behavior: select all tools
            final_selected_tools <- availableTools
        }
    } else {
        # Validate user-provided tool names
        unmatched <- setdiff(selectedTools, availableTools)
        if (length(unmatched) > 0) {
            warning("The following tools were not found and will be ignored: ",
                    paste(unmatched, collapse = ", "), call. = FALSE)
        }
        final_selected_tools <- intersect(selectedTools, availableTools)
    }
    
    if (length(final_selected_tools) == 0) {
        warning("No valid tools were selected for analysis.", call. = FALSE)
    }
    
    return(final_selected_tools)
}

#' Calculate Confusion Matrix Metrics (Internal Function)
#'
#' A lightweight, dependency-free function to calculate a suite of performance
#' metrics from prediction and reference vectors. Assumes '1' is the positive
#' class (non-coding) and '0' is the negative class (coding).
#'
#' @param predictions A numeric or integer vector of predicted labels (0 or 1).
#' @param reference A numeric or integer vector of true labels (0 or 1).
#'
#' @return A named numeric vector of performance metrics. Returns NULL if
#'   inputs are invalid.
#' @importFrom stats binom.test mcnemar.test
#' @noRd
calculateMetrics <- function(predictions, reference) {
    
    if (length(unique(predictions)) < 2) {
        warning("Prediction vector for a tool contains only one class; some metrics will be NA.",
                call. = FALSE)
    }
    
    # --- Core Confusion Matrix Components ---
    TP <- sum(predictions == 1 & reference == 1, na.rm = TRUE)
    TN <- sum(predictions == 0 & reference == 0, na.rm = TRUE)
    FP <- sum(predictions == 1 & reference == 0, na.rm = TRUE)
    FN <- sum(predictions == 0 & reference == 1, na.rm = TRUE)
    
    total_obs <- TP + TN + FP + FN
    if (total_obs == 0) return(NULL)
    
    # --- Primary Metrics ---
    Sensitivity <- TP / (TP + FN)
    Specificity <- TN / (TN + FP)
    Precision <- TP / (TP + FP)
    NegPredValue <- TN / (TN + FN)
    Accuracy <- (TP + TN) / total_obs
    F1 <- 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
    BalancedAccuracy <- (Sensitivity + Specificity) / 2
    
    # --- Statistical Metrics from caret ---
    Prevalence <- (TP + FN) / total_obs
    DetectionRate <- TP / total_obs
    DetectionPrevalence <- (TP + FP) / total_obs
    
    AccuracyNull <- max(Prevalence, 1 - Prevalence)
    
 #   accuracy_ci <- tryCatch(
 #       stats::binom.test(TP + TN, total_obs)$conf.int[1:2],
 #       error = function(e) c(NA, NA)
 #   )
    accuracy_ci <- tryCatch(
        stats::binom.test(TP + TN, total_obs)$conf.int[seq_len(2)],
        error = function(e) c(NA, NA)
    )
    AccuracyPValue <- tryCatch(
        stats::binom.test(TP + TN, total_obs, p = AccuracyNull, 
                          alternative = "greater")$p.value,
        error = function(e) NA
    )
    
    mcnemar_table <- matrix(c(TP, FN, FP, TN), nrow = 2)
    McnemarPValue <- tryCatch(
        stats::mcnemar.test(mcnemar_table)$p.value,
        error = function(e) NA
    )
    
    prob_pred_pos <- (TP + FP) / total_obs
    prob_ref_pos <- (TP + FN) / total_obs
    expected_accuracy <- (prob_pred_pos * prob_ref_pos) + 
        ((1 - prob_pred_pos) * (1 - prob_ref_pos))
    Kappa <- (Accuracy - expected_accuracy) / (1 - expected_accuracy)
    
    # --- Assemble Output Vector ---
    metrics <- c(
        "Accuracy" = Accuracy, "Kappa" = Kappa, "AccuracyLower" = accuracy_ci[1],
        "AccuracyUpper" = accuracy_ci[2], "AccuracyNull" = AccuracyNull,
        "AccuracyPValue" = AccuracyPValue, "McnemarPValue" = McnemarPValue,
        "Sensitivity" = Sensitivity, "Specificity" = Specificity,
        "Pos Pred Value" = Precision, "Neg Pred Value" = NegPredValue,
        "Precision" = Precision, "Recall" = Sensitivity, "F1" = F1,
        "Prevalence" = Prevalence, "Detection Rate" = DetectionRate,
        "Detection Prevalence" = DetectionPrevalence,
        "Balanced Accuracy" = BalancedAccuracy
    )
    
    metrics[is.nan(metrics)] <- NA
    return(metrics)
}

#' Extract a Metric Value from a CM List Object (Internal Helper)
#'
#' Retrieves a specific metric value from the simple list structure created by
#' `calculateCM`.
#'
#' @param cmObject A list containing a named numeric vector `$metrics`.
#' @param metricName The name of the metric to extract.
#' @return The numeric value of the metric, or `NA_real_` if not found.
#' @noRd
extractMetricValue <- function(cmObject, metricName) {
    value <- cmObject$metrics[metricName]
    if (is.null(value) || !is.numeric(value) || length(value) == 0) {
        NA_real_
    } else {
        value
    }
}

#' Set Colors for Plotting (Internal Helper)
#'
#' Manages color assignments for plots.
#'
#' @param colors A user-provided vector of colors, or `NULL`.
#' @param nItems The number of items that need colors.
#' @return A character vector of colors of length `nItems`.
#' @noRd
setColors <- function(colors, nItems) {
    if (is.null(colors)) {
        scales::hue_pal()(nItems)
    } else if (length(colors) < nItems) {
        warning("Fewer colors than items; recycling colors.", call. = FALSE)
        rep_len(colors, nItems)
    } else {
        colors[seq_len(nItems)]
    }
}

#' Calculate Optimal Grid Layout for Multiple Plots (Internal Helper)
#'
#' @param nPlots The number of plots to arrange.
#' @return A numeric vector `c(rows, cols)`.
#' @noRd
calculateOptimalLayout <- function(nPlots) {
    if (nPlots <= 1) return(c(1, 1))
    if (nPlots <= 2) return(c(1, 2))
    if (nPlots <= 4) return(c(2, 2))
    if (nPlots <= 6) return(c(2, 3))
    if (nPlots <= 9) return(c(3, 3))
    rows <- ceiling(sqrt(nPlots))
    cols <- ceiling(nPlots / rows)
    return(c(rows, cols))
}
