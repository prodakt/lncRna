#' Calculate Confusion Matrices and Performance Metrics for Prediction Methods
#'
#' This function calculates confusion matrices and extracts performance metrics for various
#' prediction methods compared against a reference standard. It iterates through columns
#' of prediction data, calculates confusion matrices using the `caret::confusionMatrix` function,
#' and processes/prints specified metrics based on a given threshold.
#'
#' @param bp_cmb_data data.frame. A data frame containing prediction results from different methods.
#'        Each column should represent predictions of a method and have the same number of rows
#'        as `best_pat3_data`. Column names will be used as method names.
#' @param best_pat3_data data.frame. A data frame containing the reference standard in the 'isNC' column
#'        and prediction columns that correspond to the methods in `bp_cmb_data`. Must contain a column named 'isNC'
#'        with the true class labels.
#' @param positive_class character. The class level to be considered the "positive" class when
#'        calculating confusion matrices (e.g., "1" or "TRUE"). Defaults to "1".
#' @param print_metric_threshold_methods logical. If TRUE, prints performance metric values for each method
#'        to the console, indicating if the metric meets or exceeds the specified `threshold`. Defaults to FALSE.
#' @param threshold numeric. A numeric threshold (between 0 and 1) used to determine if a method's performance
#'        on `metric_to_extract` is considered "high". Used only when `print_metric_threshold_methods = TRUE`
#'        or `return_only_high_methods = TRUE`. Defaults to 0.8.
#' @param return_only_high_methods logical. If TRUE, the function returns a list containing only the confusion matrices
#'        for methods that meet or exceed the `threshold` for the specified `metric_to_extract`. If FALSE, returns
#'        a list with confusion matrices for all methods. Defaults to FALSE.
#' @param metric_to_extract character. The name of the performance metric to extract and evaluate against the
#'        `threshold`. This metric must be a valid metric name from the output of `caret::confusionMatrix`
#'        (e.g., "Accuracy", "Sensitivity", "Specificity", "Kappa"). Defaults to "Accuracy".
#'
#' @return list.
#' If `return_only_high_methods = FALSE`, returns a named list of confusion matrix objects, where each element
#' is the confusion matrix for a method (named by method column name from `bp_cmb_data`).
#' If `return_only_high_methods = TRUE`, returns a named list containing only confusion matrices for methods
#' that have a `metric_to_extract` value greater than or equal to the specified `threshold`.
#' Returns an empty list if no methods meet the threshold, when `return_only_high_methods = TRUE`.
#'
#' @examples
#' # Assuming BP.cmb and BestPat3 data frames are already loaded
#'
#' # Example 1: Calculate confusion matrices for all columns in BP.cmb, print Accuracy, no thresholding
#' all_cms_accuracy <- calculate_cm(BP.cmb, BestPat3, print_metric_threshold_methods = TRUE, metric_to_extract = "Accuracy")
#' print(all_cms_accuracy)
#'
#' # Example 2: Calculate and return only confusion matrices with Sensitivity >= 0.75
#' high_sensitivity_cms <- calculate_cm(BP.cmb, BestPat3,
#'                                        return_only_high_methods = TRUE,
#'                                        metric_to_extract = "Sensitivity",
#'                                        threshold = 0.75,
#'                                        print_metric_threshold_methods = TRUE)
#' print(high_sensitivity_cms)
#'
#' # Example 3: Calculate with a different positive class and extract Specificity
#' specific_positive_class_cms <- calculate_cm(BP.cmb, BestPat3,
#'                                              positive_class = "0", # If '0' is your positive class
#'                                              metric_to_extract = "Specificity",
#'                                              print_metric_threshold_methods = TRUE,
#'                                              threshold = 0.9)
#' print(specific_positive_class_cms)
#'
#' # Example 4: Basic calculation, returning all confusion matrices, no printing or threshold
#' all_cms_basic <- calculate_cm(BP.cmb, BestPat3)
#' print(all_cms_basic)
#'
#' # Example 5: Using Kappa statistic as the metric to evaluate
#' kappa_metric_cms <- calculate_cm(BP.cmb, BestPat3,
#'                                   metric_to_extract = "Kappa",
#'                                   print_metric_threshold_methods = TRUE,
#'                                   threshold = 0.6,
#'                                   return_only_high_methods = TRUE)
#' print(kappa_metric_cms)
#' @import caret
#' @export
calculate_cm <- function(bp_cmb_data, best_pat3_data, positive_class = "1", print_metric_threshold_methods = FALSE, 
                         threshold = 0.8, return_only_high_methods = FALSE, metric_to_extract = "Accuracy") {
  
  # --- Input Validation ---
  validate_input_cm(bp_cmb_data, best_pat3_data, threshold, metric_to_extract)
  
  # --- Initialize lists to store results
  confusion_matrix_list <- list()
  high_metric_matrix_list <- list() # Renamed for generality
  
  # --- Get column names from bp_cmb_data (methods)
  method_columns <- colnames(bp_cmb_data)
  if (is.null(method_columns) || length(method_columns) == 0) {
    stop("Error: bp_cmb_data has no column names to copy or is empty.")
  }
  
  cat("\n --- Confusion Matrix Calculation ---\n") # header for summary output
  
  # --- Loop through each method column in bp_cmb_data
  for (selected_method in method_columns) {
    if (!selected_method %in% colnames(best_pat3_data)) {
      warning(paste("Warning: Column '", selected_method, "' not found in best_pat3_data. Skipping", sep = ""))
      next
    }
    
    # Calculate Confusion Matrix using the new function
    cm <- calculate_single_confusion_matrix(best_pat3_data, selected_method, positive_class)
    
    # Process metric, print, and update lists using the new function
    updated_lists <- process_metric_and_print(cm, selected_method, metric_to_extract, print_metric_threshold_methods, threshold, confusion_matrix_list, high_metric_matrix_list)
    confusion_matrix_list <- updated_lists$confusion_matrix_list
    high_metric_matrix_list <- updated_lists$high_metric_matrix_list
    
  }
  
  cat("\n--- End of Confusion Matrix Calculation ---\n") # Footer for summary output
  
  if (return_only_high_methods) {
    return(high_metric_matrix_list) # Return only high metric matrices if requested (renamed list)
  } else {
    return(confusion_matrix_list) # Return all matrices by default
  }
}

#######################################
# Functions Documentation (Helper Functions)

#' Validate Input Data for Confusion Matrix Calculation
#'
#' This helper function validates the input data frames and parameters
#' provided to the `calculate_cm` function, ensuring they meet the required format
#' and constraints before proceeding with confusion matrix calculations.
#'
#' @param bp_cmb_data data.frame. Input data frame containing prediction methods.
#' @param best_pat3_data data.frame. Input data frame containing reference standard ('isNC' column).
#' @param threshold numeric. Threshold value for accuracy metric comparison.
#' @param metric_to_extract character. Name of the metric to extract.
#'
#' @return No return value. Stops execution with an error message if validation fails.
#'
#' @noRd
validate_input_cm <- function(bp_cmb_data, best_pat3_data, threshold, metric_to_extract) {
  # --- Input Validation for data frames
  if(!is.data.frame(bp_cmb_data) || !is.data.frame(best_pat3_data)) {
    stop("Error: Both bp_cmb_data and best_pat3_data must be data frames.")
  }
  if(is.null(colnames(bp_cmb_data))) {
    stop("bp_cmb_data must have column names representing tools.")
  }
  if(!"isNC" %in% colnames(best_pat3_data)) {
    stop("best_pat3_data must have a column named 'isNC'.")
  }
  
  #--- Input validation for accuracy threshold
  if (!is.numeric(threshold)) {
    stop("Error: 'threshold' must be a numeric value.")
  }
  if (threshold < 0 || threshold > 1) {
    stop("Error: 'threshold' must be between 0 and 1 (for Accuracy metric comparison).")
  }
  
  # --- Input validation for metric_to_extract
  if (!is.character(metric_to_extract) || length(metric_to_extract) != 1) {
    stop("Error: 'metric_to_extract' must be a single character string (e.g., 'Accuracy', 'Kappa').")
  }
}

#' Calculate Single Confusion Matrix
#'
#' Helper function to calculate a confusion matrix for a single prediction method.
#' Uses `caret::confusionMatrix` to perform the calculation and handles potential errors.
#'
#' @param best_pat3_data data.frame. Data frame containing both prediction and reference standard.
#' @param selected_method character. Name of the column in `best_pat3_data` containing predictions for the method.
#' @param positive_class character. The class to be considered as the positive class.
#'
#' @return confusionMatrix or character.
#' Returns a confusionMatrix object if calculation is successful.
#' Returns an error message (character string) if `caret::confusionMatrix` fails.
#'
#' @noRd
calculate_single_confusion_matrix <- function(best_pat3_data, selected_method, positive_class) {
  cm <- tryCatch({ # Use tryCatch for error handling
    caret::confusionMatrix(
      data      = as.factor(best_pat3_data[, selected_method]), # Predictions from method
      reference = as.factor(best_pat3_data$isNC),          # True labels (reference)
      positive  = positive_class                           # Positive class variable
    )
  }, error = function(e) {
    # Error message if confusionMatrix fails for a method
    error_message <- paste("Error calculating confusionMatrix for method:", selected_method, "- ", e$message)
    warning(error_message) # Output a warning to the console
    return(error_message)  # Return error message
  })
  return(cm)
}

#' Process Metric, Print Results, and Update Lists
#'
#' Helper function to process the confusion matrix, extract the specified metric,
#' print results to the console (if requested), and update the lists of confusion matrices.
#'
#' @param cm confusionMatrix or character. A confusionMatrix object or an error message (character)
#'        returned from `calculate_single_confusion_matrix`.
#' @param selected_method character. Name of the method for which the confusion matrix was calculated.
#' @param metric_to_extract character. The name of the metric to extract and evaluate.
#' @param print_metric_threshold_methods logical. Flag to control printing of metric values and threshold status.
#' @param threshold numeric. Threshold value for metric comparison.
#' @param confusion_matrix_list list. List to store all confusion matrices. Modified in place.
#' @param high_metric_matrix_list list. List to store confusion matrices with high metric values. Modified in place.
#'
#' @return list. Returns a list containing the updated `confusion_matrix_list` and `high_metric_matrix_list`.
#'
#' @noRd
process_metric_and_print <- function(cm, selected_method, metric_to_extract, print_metric_threshold_methods, threshold, confusion_matrix_list, high_metric_matrix_list) {
  if (!inherits(cm, "character")) { # Check if cm is not an error message
    confusion_matrix_list[[selected_method]] <- cm # Store all CMs
    
    # --- Extract the specified metric using helper function ---
    metric_value <- extract_metric_value(cm, metric_to_extract)
    
    if (print_metric_threshold_methods) { # Check if metric threshold printing is enabled
      output_string <- paste("Method: '", selected_method, "' - ", metric_to_extract, ": ", sprintf("%.4f", metric_value), sep="") # Base output string
      
      if (!is.na(metric_value) && metric_value > threshold) { # Check if metric meets threshold (NA check added)
        output_string <- paste(output_string, " -  *** HIGH ", toupper(metric_to_extract), " (>", threshold, ") ***", sep="") # Append high metric indicator
      }
      cat(output_string, "\n")
    } else { # If print_metric_threshold_methods is FALSE, just print metric
      if (print_metric_threshold_methods) {
        cat(paste("Method: '", selected_method, "' - ", metric_to_extract, ": ", sprintf("%.4f", metric_value), "\n", sep=""))
      }
    }
    
    if (!is.na(metric_value) && metric_value >= threshold) { # Use >= for threshold comparison, NA check added
      high_metric_matrix_list[[selected_method]] <- cm # Store in separate list if high metric
    }
    
  } else {
    cat(paste("Method: '", selected_method, "' - Confusion Matrix calculation failed: ", cm, "\n", sep="")) # Indicate failure with error message
  }
  return(list(confusion_matrix_list = confusion_matrix_list, high_metric_matrix_list = high_metric_matrix_list)) # Return updated lists
}

#' Extract Metric Value from Confusion Matrix Object
#'
#' Helper function to extract a specific performance metric value from a `confusionMatrix` object.
#'
#' @param confusion_matrix_obj confusionMatrix. A confusionMatrix object from the `caret` package.
#' @param metric_name character. The name of the metric to extract (e.g., "Accuracy", "Sensitivity").
#'        Must be a valid metric name present in `confusion_matrix_obj$overall` or `confusion_matrix_obj$byClass`.
#'
#' @return numeric. The value of the extracted performance metric. Returns NA if metric is not found.
#'
#' @noRd
extract_metric_value <- function(confusion_matrix_obj, metric_name) {
  metric_value <- NA # Default value if metric not found
  
  if (metric_name %in% names(confusion_matrix_obj$overall)) {
    metric_value <- confusion_matrix_obj$overall[metric_name]
  } else if (metric_name %in% names(confusion_matrix_obj$byClass)) {
    metric_value <- confusion_matrix_obj$byClass[metric_name]
  } else {
    stop(paste("Warning: Metric '", metric_name, "' not found in confusion matrix. Plese enter correct metric.", sep=""))
  }
  return(metric_value)
}

