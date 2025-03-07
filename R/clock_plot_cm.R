#' Create Clock Plots for Confusion Matrix Metrics
#'
#' Generates clock plots (also known as circular bar plots) to visualize performance metrics
#' derived from confusion matrices. This function supports both single and multiple plot layouts,
#' allowing for comparison of different methods or a detailed view of individual method performance.
#'
#' @param cm_list list. A named list of confusion matrix objects, typically created using `caret::confusionMatrix`.
#'        Names of the list elements should correspond to the method names.
#' @param methods character vector, optional. A vector of method names to include in the plot.
#'        If NULL, the function will prompt the user to interactively select methods from `cm_list`.
#' @param metrics character vector. A vector of metric names to be displayed in the clock plot.
#'        Defaults to c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall").
#'        Valid metrics include "Accuracy", "Sensitivity", "Specificity", "Precision", "Recall",
#'        "Pos Pred Value", "Neg Pred Value", "F1", "Prevalence", "Detection Rate",
#'        "Detection Prevalence", "Balanced Accuracy", and any other metrics available in the
#'        `confusionMatrix` object's `overall` or `byClass` slots.
#' @param plot_title character, optional. The title of the clock plot. Defaults to "Clock Plot of Metrics".
#' @param colors character vector, optional. A vector of colors to use for each method in the plot.
#'        If NULL, a default color palette from `scales::hue_pal()` is used. If fewer colors are provided
#'        than methods, colors will be recycled.
#' @param layout character, optional. Specifies the layout of the plots.
#'        "single": Generates a single clock plot comparing all selected methods.
#'        "multiple": Generates a series of clock plots, one for each selected method.
#'        Defaults to "single".
#' @param ... Additional arguments to be passed to the underlying plotting functions (ggplot2).
#'
#' @return  ggplot object or patchwork object.
#' Returns a ggplot object when `layout = "single"`.
#' Returns a patchwork object (arranging multiple ggplot plots) when `layout = "multiple"`.
#' These objects can be further customized or printed.
#'
#' @examples
#' # Assuming 'all_cms' is a named list of confusionMatrix objects
#'
#' # Example 1: Single clock plot with interactive method selection
#' clock_plot_cm(all_cms, plot_title = "Single Clock Plot Comparison")
#'
#' # Example 2: Multiple clock plots, one for each method
#' clock_plot_cm(all_cms, layout = "multiple", plot_title = "Method-Specific Clock Plots")
#'
#' # Example 3: Single clock plot for specific methods and metrics
#' clock_plot_cm(all_cms, methods = c("CPC2+CPAT", "TargetScan"),
#'                metrics = c("Accuracy", "Sensitivity", "Specificity"),
#'                plot_title = "Custom Single Clock Plot")
#'
#' # Example 4: Multiple clock plots with custom colors
#' clock_plot_cm(all_cms, methods = names(all_cms)[1:3], layout = "multiple",
#'                plot_title = "Clock Plots with Custom Colors",
#'                colors = c("red", "blue", "green"))
#'
#' # Example 5: Single plot with different metrics
#' clock_plot_cm(all_cms, metrics = c("F1", "Balanced Accuracy", "Precision", "Recall"),
#'                plot_title = "Clock Plot - Alternative Metrics")
#'
#' @export
clock_plot_cm <- function(cm_list, methods = NULL, 
                          metrics = c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall"),
                          plot_title = "Clock Plot of Metrics", colors = NULL, layout = "single",...) {
  layout <- match.arg(layout, choices = c("single", "multiple"))
  
  # Validate and prepare inputs
  validated_inputs <- validate_inputs_clock_plot(cm_list, methods, metrics)
  cm_list <- validated_inputs$cm_list
  methods <- validated_inputs$methods
  metrics <- validated_inputs$metrics
  
  # Prepare clock plot data
  clock_data <- prepare_clock_data(cm_list, methods, metrics)
  
  # Set colors
  colors <- set_colors(colors, length(methods))
  
  # Draw clock plot based on layout
  plot_clock(clock_data, plot_title, colors, layout, ...)
  
}

#######################################
# Helper Functions Documentation

#' Validate Inputs for Clock Plot Functions
#'
#' Validates the input list of confusion matrices, selected methods, and metrics
#' to ensure they are in the correct format for clock plot generation. When methods
#' are not provided, it offers an interactive prompt with "All methods" as the first option.
#'
#' @param cm_list list. A named list of confusion matrix objects.
#' @param methods character vector, optional. Vector of method names to be plotted, or NULL for interactive selection.
#' @param metrics character vector. Vector of metric names to be plotted.
#'
#' @return list. Returns a list containing validated `cm_list`, `methods`, and `metrics`.
#'
#' @noRd
validate_inputs_clock_plot <- function(cm_list, methods, metrics) {
  # Check if cm_list is a valid named list
  if (!is.list(cm_list) || is.null(names(cm_list)) || length(names(cm_list)) == 0) {
    stop("cm_list must be a named list of confusion matrix objects.")
  }
  available_methods <- names(cm_list)
  
  # Define valid metrics (assuming these are typical confusion matrix metrics)
  valid_metrics <- unique(c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall",
                            "Pos Pred Value", "Neg Pred Value", "F1", "Prevalence",
                            "Detection Rate", "Detection Prevalence", "Balanced Accuracy",
                            names(cm_list[[1]]$overall), names(cm_list[[1]]$byClass)))
  
  # Interactive method selection if methods is NULL
  if (is.null(methods)) {
    cat("Available methods in cm_list:\n")
    cat("1. All methods\n")
    for (i in seq_along(available_methods)) {
      cat(paste0(i + 1, ". ", available_methods[i], "\n"))
    }
    selected_indices_str <- readline(prompt = "Select method numbers (comma-separated, e.g., 1 or 2,3): ")
    selected_indices <- as.numeric(strsplit(selected_indices_str, ",")[[1]])
    
    # Validate input
    if (any(is.na(selected_indices))) {
      stop("Error: Invalid input. Please enter numbers separated by commas.")
    }
    if (any(selected_indices < 1 | selected_indices > length(available_methods) + 1)) {
      stop(sprintf("Error: Method index out of range. Please select indices between 1 and %d.", 
                   length(available_methods) + 1))
    }
    
    # Handle "All methods" or specific selections
    if (1 %in% selected_indices) {
      methods <- available_methods
      cat("You selected: All methods\n")
    } else {
      # Adjust indices since "All methods" is 1, and individual methods start from 2
      adjusted_indices <- selected_indices - 1
      methods <- available_methods[adjusted_indices]
      cat(sprintf("You selected methods: '%s'\n", paste(methods, collapse = "', '")))
    }
  } else if (!all(methods %in% available_methods)) {
    # Validate explicitly provided methods
    stop(sprintf("Methods '%s' not in cm_list. Available: '%s'.",
                 paste(methods[!methods %in% available_methods], collapse = "', '"),
                 paste(available_methods, collapse = "', '")))
  }
  
  # Validate metrics
  if (!all(metrics %in% valid_metrics)) {
    stop(sprintf("Invalid metrics '%s'. Valid options: '%s'.",
                 paste(metrics[!metrics %in% valid_metrics], collapse = "', '"),
                 paste(valid_metrics, collapse = "', '")))
  }
  
  # Return validated inputs
  list(cm_list = cm_list, methods = methods, metrics = metrics)
}

#' Prepare Data for Clock Plot
#'
#' Transforms a list of confusion matrix objects into a data frame suitable for creating
#' clock plots with ggplot2. It extracts the specified metrics for each method and structures
#' the data for plotting.
#'
#' @param cm_list list. A named list of confusion matrix objects.
#' @param methods character vector. Vector of method names to include in the data.
#' @param metrics character vector. Vector of metric names to extract and include in the data.
#'
#' @return data.frame. A data frame structured for clock plot generation, with columns for Method, Metric, and Value.
#'
#' @noRd
prepare_clock_data <- function(cm_list, methods, metrics) {
  clock_data_list <- lapply(methods, function(method) {
    values <- sapply(metrics, function(metric) {
      extract_metric_value(cm_list[[method]], metric)
    })
    data.frame(
      Method = rep(method, length(metrics)),
      Metric = metrics,
      Value = values,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, clock_data_list)
}

#' Set Colors for Clock Plots
#'
#' Manages color assignment for clock plots. If colors are not provided, it uses a default
#' palette. If fewer colors are provided than methods, it recycles colors with a warning.
#'
#' @param colors character vector, optional. User-provided color vector, or NULL to use default colors.
#' @param n_methods integer. The number of methods for which colors are needed.
#'
#' @return character vector. A vector of colors of length `n_methods`.
#'
#' @noRd
set_colors <- function(colors, n_methods) {
  if (is.null(colors)) {
    scales::hue_pal()(n_methods)
  } else if (length(colors) < n_methods) {
    warning("Fewer colors than methods; recycling colors.")
    rep_len(colors, n_methods)
  } else {
    colors[1:n_methods]
  }
}

#' Dispatch Clock Plot Drawing Based on Layout
#'
#' Determines whether to draw a single or multiple clock plots based on the specified layout.
#' This function acts as a dispatcher to call the appropriate plotting function.
#'
#' @param clock_data data.frame. Data frame prepared by `prepare_clock_data`.
#' @param plot_title character. The title of the plot(s).
#' @param colors character vector. Vector of colors for methods.
#' @param layout character. Specifies "single" or "multiple" layout.
#' @param ... Additional arguments passed to `draw_single_clock_plot` or `draw_multiple_clock_plots`.
#'
#' @return ggplot object or list of patchwork objects. Returns the plot object(s) generated by the called function.
#'
#' @noRd
plot_clock <- function(clock_data, plot_title, colors, layout, ...) {
  if (layout == "single") {
    draw_single_clock_plot(clock_data, plot_title, colors, ...)
  } else {
    draw_multiple_clock_plots(clock_data, plot_title, colors, ...)
  }
}

#' Draw a Single Clock Plot
#'
#' Generates a single clock plot comparing performance metrics across different methods using ggplot2.
#'
#' @param clock_data data.frame. Data frame prepared by `prepare_clock_data`.
#' @param plot_title character. The title of the plot.
#' @param colors character vector. Vector of colors for methods.
#' @param ... Additional arguments passed to `ggplot2` geom or theme functions.
#'
#' @return ggplot object. Returns a ggplot object representing the single clock plot.
#'
#' @noRd
draw_single_clock_plot <- function(clock_data, plot_title, colors, ...) {
  # Base plot
  p <- ggplot(clock_data, aes(x = Metric, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
    coord_polar(start = 0) +
    ylim(0, 1.1) +  # Extend y-axis slightly for labels
    ggtitle(plot_title) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    geom_text(aes(y = Value + 0.03, label = sprintf("%.3f", Value)),
              position = position_dodge(width = 0.9), vjust = 0, size = 3)
  
  print(p)
}

#' Draw Multiple Clock Plots (One per Method)
#'
#' Generates a series of clock plots with performance metrics for each method, ensuring metric names are fully visible.
#'
#' @param clock_data data.frame. Data frame prepared by `prepare_clock_data`.
#' @param plot_title character. The base title for the plots, displayed on the first line.
#' @param colors character vector. Vector of colors for methods.
#' @param ... Additional arguments passed to `ggplot2` geom or theme functions.
#'
#' @return list of patchwork objects. Each object is a grid of up to 4 clock plots with two-line titles.
#'
#' @noRd
draw_multiple_clock_plots <- function(clock_data, plot_title, colors, ...) {
  library(stringr)  # For str_wrap
  methods <- unique(clock_data$Method)
  num_methods <- length(methods)
  if (num_methods <= 0) stop("No methods available to plot.")
  
  # Reduce to 2x2 grid for more space per plot
  max_plots_per_grid <- 6
  n_grids <- ceiling(num_methods / max_plots_per_grid)
  
  # Split methods into chunks of up to 4
  method_chunks <- split(methods, rep(1:n_grids, each = max_plots_per_grid, length.out = num_methods))
  
  # Create a list to hold all grid plots
  all_plots <- list()
  
  for (grid in 1:n_grids) {
    current_methods <- method_chunks[[grid]]
    
    # Create plots for the current grid
    plot_list <- lapply(current_methods, function(method) {
      method_data <- clock_data[clock_data$Method == method, ]
      # Wrap long metric names
      method_data$Metric <- str_wrap(method_data$Metric, width = 10)
      title_text <- paste(plot_title, "\n", method)
      ggplot(method_data, aes(x = Metric, y = Value)) +
        geom_bar(stat = "identity", fill = alpha(colors[which(methods == method)], 0.5)) +
        coord_polar(start = 0) +
        ylim(0, 1.3) +  # Extend radial limit further
        ggtitle(title_text) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
          axis.text.y = element_blank(),
          # Remove default x-axis labels
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
          plot.margin = unit(c(1, 1, 2, 1), "cm")  # Extra bottom margin
        ) +
        # Add custom labels outside the plot
        geom_text(aes(label = Metric, x = Metric, y = 1.25), 
                  angle = 0, hjust = 0.5, vjust = 0, size = 3.5) +
        geom_text(aes(y = Value + 0.05, label = sprintf("%.4f", Value)),
                  position = position_dodge(width = 0.9), vjust = 0, size = 3)
    })
    
    # Arrange plots in a 2x2 grid
    grid_layout <- patchwork::wrap_plots(plot_list) +
      plot_layout(guides = "collect") & 
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    all_plots[[grid]] <- grid_layout
  }
  
  return(all_plots)
}

#' Extract Metric Value from Confusion Matrix Object
#'
#' Helper function to extract a specific performance metric value from a `confusionMatrix` object.http://127.0.0.1:30639/graphics/plot_zoom_png?width=1920&height=1017
#'  Reused from radar plot functions for consistency.
#'
#' @param confusion_matrix_obj confusionMatrix. A confusionMatrix object from the `caret` package.
#' @param metric_name character. The name of the metric to extract (e.g., "Accuracy", "Sensitivity").
#'        Must be a valid metric name present in `confusion_matrix_obj$overall` or `confusion_matrix_obj$byClass`.
#'
#' @return numeric. The value of the extracted performance metric.
#'
#' @noRd
extract_metric_value <- function(confusion_matrix_obj, metric_name) {
  if (metric_name %in% names(confusion_matrix_obj$overall)) {
    confusion_matrix_obj$overall[metric_name]
  } else if (metric_name %in% names(confusion_matrix_obj$byClass)) {
    confusion_matrix_obj$byClass[metric_name]
  } else {
    stop(sprintf("Metric  not found in confusion matrix.", metric_name))
  }
}

