#' Create Radar Plots for Confusion Matrix Metrics
#'
#' Generates radar plots (also known as spider or star plots) to visualize and compare
#' performance metrics from confusion matrices of different methods. Supports single plot
#' to compare multiple methods or multiple plots to evaluate each method individually.
#'
#' @param cm_list list. A named list of confusion matrix objects, typically output from
#'        `caret::confusionMatrix`. List names should correspond to the method names.
#' @param methods character vector, optional. Vector of method names to be included in the plot.
#'        If NULL, the function will provide an interactive prompt to select methods from `cm_list`.
#' @param metrics character vector. Metrics to be visualized on the radar plot axes.
#'        Defaults to c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall").
#'        Valid options include metrics from `confusionMatrix` output like "Accuracy", "Sensitivity",
#'        "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1",
#'        "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy",
#'        and any other metric available in `confusion_matrix_obj$overall` or `$byClass`.
#' @param plot_title character, optional. Title of the radar plot. Defaults to "Radar Plot of Metrics".
#' @param colors character vector, optional. Colors to be used for each method's radar line and fill.
#'        If NULL, default colors from `scales::hue_pal()` are used. If fewer colors than methods are provided,
#'        colors are recycled.
#' @param layout character, optional. Specifies the plot layout:
#'        "single": Creates a single radar plot comparing all selected methods.
#'        "multiple": Generates a series of radar plots, one for each method.
#'        Defaults to "single".
#' @param display_area logical, optional. If TRUE, displays the polygon area for each method,
#'        either in the legend for single plots or below each plot in multiple layouts. Defaults to FALSE.
#' @param display_fill logical, optional. If TRUE, fills the radar chart polygons with a semi-transparent color
#'        corresponding to each method. If FALSE, polygons are not filled. Defaults to TRUE.
#' @param save_data logical, optional. If TRUE, saves the radar data frame to a CSV file. Defaults to FALSE.
#' @param file_name character, optional. The name of the file to save the radar data frame. Required if `save_data` is TRUE.
#' @param ... Additional arguments passed to the `fmsb::radarchart` function for customization
#'        of the radar plot appearance.
#'
#' @return ggplot object or NULL.
#' Returns `NULL` as the function primarily produces plots as side effects.
#' When `layout = "single"`, the plot is directly printed.
#' When `layout = "multiple"`, multiple plots are generated and arranged using `par(mfrow=...)`.
#'
#' @examples
#' # Assuming 'all_cms' is a named list of confusionMatrix objects
#'
#' # Example 1: Single radar plot comparing methods, area in legend, no fill
#' radar_plot_cm(all_cms, plot_title = "Comparison Radar Plot", display_area = TRUE, display_fill = FALSE)
#'
#' # Example 2: Multiple radar plots, one per method, with area and fill
#' radar_plot_cm(all_cms, layout = "multiple", plot_title = "Radar Plots per Method", display_area = TRUE, display_fill = TRUE)
#'
#' # Example 3: Single plot, specific methods and metrics, custom colors
#' radar_plot_cm(all_cms, methods = c("Method1", "Method3"),
#'                metrics = c("Accuracy", "F1", "Specificity"),
#'                plot_title = "Custom Metric Radar Plot",
#'                colors = c("skyblue", "salmon"))
#'
#' # Example 4: Multiple plots, different metrics, no area display
#' radar_plot_cm(all_cms, layout = "multiple",
#'                metrics = c("Precision", "Recall", "Balanced Accuracy"),
#'                plot_title = "Multiple Radar Plots - Alt Metrics",
#'                display_area = FALSE)
#'
#' # Example 5: Single plot, interactive method selection will be prompted if methods=NULL
#' # radar_plot_cm(all_cms, plot_title = "Interactive Method Radar Plot") # Uncomment to run interactively
#'
#' @import fmsb
#' @import scales
#' @export
radar_plot_cm <- function(cm_list, methods = NULL, metrics = c("Accuracy", "Sensitivity", "Specificity","Precision","Recall"),
                          plot_title = "Radar Plot of Metrics", colors = NULL, layout = "single",
                          display_area = FALSE, display_fill = T, save_data = F, file_name = NULL, ...) {
  
  layout <- match.arg(layout, choices = c("single", "multiple"))
  
  # Validate and prepare inputs
  validated_inputs <- validate_inputs_radar_plot(cm_list, methods, metrics)
  cm_list <- validated_inputs$cm_list
  methods <- validated_inputs$methods
  metrics <- validated_inputs$metrics
  
  # Prepare radar plot data
  radar_data <- prepare_radar_data(cm_list, methods, metrics)
  
  # Save radar data if requested
  if (save_data) {
    if (is.null(file_name) || !is.character(file_name)) {
      stop("Please provide a valid file name when save_data is TRUE.")
    }
    write.csv2(radar_data, file = file_name, row.names = TRUE)
    message("Radar data frame saved to ", file_name)
  }
  
  # Set colors
  colors <- set_colors(colors, length(methods))
  
  # Draw radar plot based on layout
  plot_radar(radar_data, plot_title, colors, layout, display_area, display_fill, ...)
}

#######################################
# Helper Functions Documentation

#' Validate Inputs for Radar Plot Functions
#'
#' Validates the inputs for radar plot functions, ensuring that `cm_list` is a properly
#' formatted list of confusion matrices, methods are valid, and metrics are recognized.
#' Handles interactive method selection if `methods` is NULL, with "All methods" as the first option.
#'
#' @param cm_list list. Named list of confusion matrix objects.
#' @param methods character vector, optional. Method names to be plotted, or NULL for interactive selection.
#' @param metrics character vector. Metric names to be plotted.
#'
#' @return list. Returns a list containing validated `cm_list`, `methods`, and `metrics`.
#'
#' @noRd
validate_inputs_radar_plot <- function(cm_list, methods, metrics) {
  # Check if cm_list is a valid named list
  if (!is.list(cm_list) || is.null(names(cm_list)) || length(names(cm_list)) == 0) {
    stop("cm_list must be a named list of confusion matrix objects.")
  }
  available_methods <- names(cm_list)
  
  # Define valid metrics based on standard names and confusion matrix fields
  valid_metrics <- unique(c("Accuracy", "Sensitivity", "Specificity", "Pos Pred Value",
                            "Neg Pred Value", "Precision", "Recall", "F1", "Prevalence",
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
    
    # Handle selections
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
    # Validate provided methods if not NULL
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

#' Prepare Data for Radar Plot
#'
#' Creates a data frame in the format required by `fmsb::radarchart` from a list of
#' confusion matrix objects, selected methods, and metrics. Scales metric values to [0, 1] if necessary.
#'
#' @param cm_list list. Named list of confusion matrix objects.
#' @param methods character vector. Method names to be included in the radar plot data.
#' @param metrics character vector. Metrics to be included in the radar plot data.
#'
#' @return data.frame. A data frame ready for use with `fmsb::radarchart`.
#'
#' @noRd
prepare_radar_data <- function(cm_list, methods, metrics) {
  radar_df <- data.frame(row.names = c("Max", "Min", methods))
  for (metric in metrics) {
    values <- sapply(methods, function(m) extract_metric_value(cm_list[[m]], metric))
    if (any(values < 0 | values > 1, na.rm = TRUE)) {
      warning(sprintf("Metric '%s' has values outside [0,1]; scaling to fit.", metric))
      values <- scales::rescale(values, to = c(0, 1))
    }
    radar_df[[metric]] <- c(1, 0, values)  # Max=1, Min=0, then method values
  }
  radar_df
}

#' Set Colors for Radar Plots
#'
#' Manages color assignments for radar plots, using default palette if no colors are provided,
#' and recycling colors if fewer are provided than methods.
#'
#' @param colors character vector, optional. User-provided colors, or NULL for default.
#' @param n_methods integer. Number of methods to color.
#'
#' @return character vector. Color vector of length `n_methods`.
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

#' Dispatch Radar Plot Drawing Based on Layout
#'
#' Determines whether to draw a single combined radar plot or multiple individual plots
#' based on the specified layout. Acts as a dispatcher to call the appropriate drawing function.
#'
#' @param radar_data data.frame. Data frame prepared by `prepare_radar_data`.
#' @param plot_title character. Title of the plot(s).
#' @param colors character vector. Colors for methods.
#' @param layout character. "single" or "multiple" layout.
#' @param display_area logical. Whether to display polygon area.
#' @param display_fill logical. Whether to fill radar polygons.
#' @param ... Additional arguments passed to `draw_single_radar_plot` or `draw_multiple_radar_plots`.
#'
#' @return NULL. Plotting is a side effect; function returns NULL.
#'
#' @noRd
plot_radar <- function(radar_data, plot_title, colors, layout, display_area, display_fill, ...) {
  if (layout == "single") {
    draw_single_radar_plot(radar_data, plot_title, colors, display_area, display_fill, ...)
  } else {
    draw_multiple_radar_plots(radar_data, plot_title, colors, display_area, display_fill, ...)
  }
}

#' Draw Single Radar Plot
#'
#' Generates a single radar plot comparing metrics across different methods using `fmsb::radarchart`.
#' Supports display of polygon area in the legend and polygon fill.
#'
#' @param radar_data data.frame. Data frame prepared by `prepare_radar_data`.
#' @param plot_title character. Title of the plot.
#' @param colors character vector. Colors for methods.
#' @param display_area logical. Whether to display polygon area in the legend.
#' @param display_fill logical. Whether to fill radar polygons.
#' @param ... Additional arguments passed to `fmsb::radarchart`.
#'
#' @return NULL. Plotting is a side effect; function returns NULL.
#'
#' @noRd
draw_single_radar_plot <- function(radar_data, plot_title, colors, display_area, display_fill, ...) {
  if (display_fill) {
    plot_fill <- scales::alpha(colors, 0.5)
  } else {
    plot_fill <- NA
  }
  fmsb::radarchart(
    radar_data,
    pfcol = plot_fill,  # Semi-transparent fill
    pcol = colors,                       # Line colors
    plwd = 2,                            # Line width
    plty = 1,                            # Line type
    cglcol = "grey",                     # Grid color
    cglty = 1,                           # Grid line type
    axislabcol = "black",                 # Axis label color
    caxislabels = seq(0, 1, 0.2),        # Axis labels
    cglwd = 0.8,                         # Grid line width
    axistype = 1,                        # Axis type
    seg = 4,                             # Number of segments
    title = plot_title,
    ...
  )
  
  # Prepare legend
  legend_text <- rownames(radar_data)[-c(1, 2)]  # Exclude Max and Min
  if (display_area) {
    areas <- sapply(rownames(radar_data)[-(1:2)], function(method) {
      area_val <- calculate_radar_polygon_area(radar_data, method, colnames(radar_data))
      paste0(method, " (Area=", sprintf("%.4f", area_val), ")")
    })
    legend_text <- areas
  }
  
  # Position legend inside plot area
  legend("bottom", legend = legend_text, bty = "n", pch = 20, col = colors,
         text.col = "black",  cex = 0.8,      # Reduced legend text size
         pt.cex = 1.5, xpd = T)   # Slightly reduced point size in legend
}

#' Draw Multiple Radar Plots (One per Method)
#'
#' Generates multiple radar plots, each displaying metrics for a single method. Arranges plots in a grid layout,
#' splitting across multiple pages if the number of methods exceeds the maximum plots per page.
#'
#' @param radar_data data.frame. Data frame prepared by `prepare_radar_data`.
#' @param plot_title character. Base title for plots; method name is appended to each.
#' @param colors character vector. Colors for methods.
#' @param display_area logical. Whether to display polygon area below each plot.
#' @param display_fill logical. Whether to fill radar polygons.
#' @param ... Additional arguments passed to `fmsb::radarchart`.
#'
#' @return NULL. Plotting is a side effect; function returns NULL.
#'
#' @noRd
draw_multiple_radar_plots <- function(radar_data, plot_title, colors, display_area, display_fill, ...) {
  num_methods <- nrow(radar_data) - 2  # Number of methods is rows minus Max/Min
  if (num_methods <= 0) stop("No methods available to plot.")
  
  max_plots_per_page <- 9  # Limit to 9 plots per page (e.g., 3x3 grid)
  n_pages <- ceiling(num_methods / max_plots_per_page)  # Calculate total pages needed
  
  methods_to_plot <- rownames(radar_data)[3:nrow(radar_data)]  # Method names
  start_idx <- 1
  
  # Save original graphical parameters
  opar <- par(no.readonly = TRUE)
  
  # Loop through pages
  for (page in 1:n_pages) {
    # Determine the number of plots for this page
    end_idx <- min(start_idx + max_plots_per_page - 1, num_methods)
    num_plots_this_page <- end_idx - start_idx + 1
    current_methods <- methods_to_plot[start_idx:end_idx]
    
    # Set layout for this page with smaller margins
    layout_dims <- calculate_optimal_layout(num_plots_this_page)
    par(mar = c(3, 3, 3, 1), mfrow = layout_dims)  # Reduced margins
    
    # Plot each method in this chunk
    for (method_name in current_methods) {
      method_idx <- which(rownames(radar_data) == method_name)
      method_data <- radar_data[c(1, 2, method_idx), ]
      if (display_fill) {
        plot_fill <- scales::alpha(colors[method_idx - 2], 0.5)
      } else {
        plot_fill <- NA
      }
      
      fmsb::radarchart(
        method_data,
        pfcol = plot_fill,
        pcol = colors[method_idx - 2],
        plwd = 2,
        plty = 1,
        cglcol = "grey",
        cglty = 1,
        axislabcol = "black",
        caxislabels = seq(0, 1, 0.2),
        cglwd = 0.8,
        axistype = 1,
        seg = 4,
        title = paste(plot_title, "- Method:", method_name),
        ...
      )
      
      if (display_area) {
        polygon_area <- calculate_radar_polygon_area(radar_data, method_name, colnames(radar_data))
        mtext(paste("Area =", sprintf("%.4f", polygon_area)), side = 1, line = 1,
              cex = 0.8, col = colors[method_idx - 2])
      }
    }
    
    # Move to the next chunk of methods
    start_idx <- end_idx + 1
  }
  
  # Restore original graphical parameters
  par(opar)
}

#' Calculate Radar Polygon Area
#'
#' Calculates the area of the polygon formed by the radar chart for a given method.
#' Used to quantify the overall performance represented by the radar plot.
#'
#' @param radar_data data.frame. Data frame prepared by `prepare_radar_data`.
#' @param method_name character. Name of the method for which to calculate the polygon area.
#' @param metrics character vector. Metrics used to create the radar plot.
#'
#' @return numeric. The calculated area of the radar plot polygon.
#'
#' @noRd
calculate_radar_polygon_area <- function(radar_data, method_name, metrics) {
  method_values <- as.numeric(radar_data[rownames(radar_data) == method_name, metrics])
  n_metrics <- length(method_values)
  angles <- seq(0, 2 * pi, length.out = n_metrics + 1)[1:n_metrics]
  x_coords <- method_values * cos(angles)
  y_coords <- method_values * sin(angles)
  
  n <- length(x_coords)
  area <- 0.5 * abs(sum(x_coords[1:n] * y_coords[c(2:n, 1)]) - sum(y_coords[1:n] * x_coords[c(2:n, 1)]))
  area
}

#' Calculate Optimal Layout for Multiple Plots
#'
#' Determines the optimal grid layout (rows and columns) for arranging a given number of plots.
#' This function is reused from clock plot functions for consistent layout management.
#'
#' @param n_plots integer. The number of plots to arrange.
#'
#' @return integer vector. A vector of length 2, specifying the number of rows and columns (c(rows, cols)).
#'
#' @noRd
calculate_optimal_layout <- function(n_plots) {
  if (n_plots <= 2) return(c(1, n_plots))
  if (n_plots <= 4) return(c(2, 2))
  if (n_plots <= 6) return(c(2, 3))
  rows <- ceiling(sqrt(n_plots))
  cols <- ceiling(n_plots / rows)
  c(rows, cols)
}

#' Extract Metric Value from Confusion Matrix Object
#'
#' Helper function to extract a specific performance metric value from a `confusionMatrix` object.
#' Reused across plotting functions for consistency in metric extraction.
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
    stop(sprintf("Metric '%s' not found in confusion matrix.", metric_name))
  }
}
