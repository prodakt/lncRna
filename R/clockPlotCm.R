#' Create Clock Plots for Confusion Matrix Metrics
#'
#' Generates circular bar plots to visualize performance metrics derived from
#' confusion matrices. Requires the 'ggplot2', 'patchwork', and 'stringr'
#' packages.
#'
#' @param cmList A named list of `confusionMatrix` objects. Names of the list
#'   elements correspond to the method names.
#' @param methods An optional character vector of method names to include in the
#'   plot. If `NULL` (default) and the session is interactive, the user will be
#'   prompted to select methods. In a non-interactive session, `NULL` results
#'   in all methods being plotted.
#' @param metrics A character vector of metric names to display. Defaults to
#'   c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall").
#' @param plotTitle An optional title for the clock plot.
#' @param colors An optional character vector of colors for each method. If
#'   `NULL`, a default color palette is used.
#' @param layout Specifies the layout: "single" for a combined plot or
#'   "multiple" for one plot per method. Defaults to "single".
#'
#' @return A `ggplot` or `patchwork` object, or `invisible(NULL)` if required
#'   packages are not installed.
#'
#' @export
#' @examples
#' if (requireNamespace("caret", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE) &&
#'     requireNamespace("stringr", quietly = TRUE)) {
#'
#'   # --- 1. Create a sample list of confusionMatrix objects ---
#'   set.seed(42)
#'   ref <- factor(sample(c("A", "B"), 100, replace = TRUE))
#'   pred1 <- factor(sample(c("A", "B"), 100, replace = TRUE))
#'   pred2 <- factor(sample(c("A", "B"), 100, replace = TRUE))
#'   sampleCmList <- list("MethodA" = caret::confusionMatrix(pred1, ref),
#'                        "MethodB" = caret::confusionMatrix(pred2, ref))
#'
#'   # --- 2. Run non-interactive examples ---
#'   p1 <- clockPlotCm(
#'     cmList = sampleCmList,
#'     methods = "MethodA",
#'     plotTitle = "Performance of MethodA"
#'   )
#'   # To display the plot in an interactive session, simply type `p1`
#'
#' } else {
#'   message("Skipping clock plot examples because dependencies are not installed.")
#' }
clockPlotCm <- function(cmList, methods = NULL,
                        metrics = c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall"),
                        plotTitle = "Clock Plot of Metrics", colors = NULL, layout = "single") {

  # --- Check for suggested packages ---
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE) ||
      !requireNamespace("stringr", quietly = TRUE)) {
    message("Packages 'ggplot2', 'patchwork', and 'stringr' are required for clockPlotCm().")
    message("Please install them using: BiocManager::install(c('ggplot2', 'patchwork', 'stringr'))")
    return(invisible(NULL))
  }

  layout <- match.arg(layout, choices = c("single", "multiple"))
  validateClockPlotInputs(cmList, metrics)
  availableMethods <- names(cmList)

  if (is.null(methods)) {
    selectedMethods <- selectMethodsInteractively(availableMethods)
  } else {
    stopifnot("Some 'methods' are not in 'cmList'" = all(methods %in% availableMethods))
    selectedMethods <- methods
  }
  if (is.null(selectedMethods) || length(selectedMethods) == 0) {
    warning("No valid methods selected. Returning NULL.")
    return(invisible(NULL))
  }

  clockData <- prepareClockData(cmList, selectedMethods, metrics)
  plotColors <- setColors(colors, length(selectedMethods))
  dispatchClockPlot(clockData, plotTitle, plotColors, layout)
}

# --- Helper Functions ---

#' Interactively Select Methods or Default to All
#' @param availableMethods A character vector of method names.
#' @return A character vector of selected method names.
#' @noRd
selectMethodsInteractively <- function(availableMethods) {
  if (!interactive()) {
    return(availableMethods)
  }

  message("Available methods:")
  choices <- c("All methods", availableMethods)
  choicesText <- paste0(seq_along(choices), ". ", choices)
  message(paste(choicesText, collapse = "\n"))

  prompt <- "Select method numbers (comma-separated, e.g., 1 or 2,3): "
  userInput <- readline(prompt = prompt)
  selectedIndices <- suppressWarnings(as.numeric(strsplit(userInput, ",")[[1]]))

  if (any(is.na(selectedIndices)) || length(selectedIndices) == 0) {
    warning("Invalid input. Defaulting to all methods.")
    return(availableMethods)
  }
  if (1 %in% selectedIndices) {
    return(availableMethods)
  }

  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(choices)]
  if (length(validIndices) == 0) {
    warning("No valid methods selected. Returning NULL.")
    return(NULL)
  }
  return(availableMethods[validIndices - 1])
}

#' Validate Inputs for Clock Plot Functions
#' @noRd
validateClockPlotInputs <- function(cmList, metrics) {
  if (!is.list(cmList) || is.null(names(cmList)) || length(names(cmList)) == 0) {
    stop("cmList must be a named list of confusion matrix objects.")
  }
  validMetrics <- unique(c("Accuracy", "Sensitivity", "Specificity", "Precision",
                           "Recall", "Pos Pred Value", "Neg Pred Value", "F1", "Prevalence",
                           "Detection Rate", "Detection Prevalence", "Balanced Accuracy",
                           names(cmList[[1]]$overall), names(cmList[[1]]$byClass)))
  if (!all(metrics %in% validMetrics)) {
    stop(sprintf("Invalid metrics '%s'.",
                 paste(metrics[!metrics %in% validMetrics], collapse = "', '")))
  }
}

#' Prepare Data for Clock Plot
#' @noRd
prepareClockData <- function(cmList, methods, metrics) {
  clockDataList <- lapply(methods, function(method) {
    values <- sapply(metrics, function(metric) {
      extractMetricValue(cmList[[method]], metric)
    })
    data.frame(
      Method = method,
      Metric = metrics,
      Value = as.numeric(values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, clockDataList)
}

#' Set Colors for Clock Plots
#' @noRd
setColors <- function(colors, nMethods) {
  if (is.null(colors)) {
    scales::hue_pal()(nMethods)
  } else if (length(colors) < nMethods) {
    warning("Fewer colors than methods; recycling colors.")
    rep_len(colors, nMethods)
  } else {
    colors[seq_len(nMethods)]
  }
}

#' Dispatch Clock Plot Drawing Based on Layout
#' @noRd
dispatchClockPlot <- function(clockData, plotTitle, colors, layout) {
  if (layout == "single") {
    drawSingleClockPlot(clockData, plotTitle, colors)
  } else {
    drawMultipleClockPlots(clockData, plotTitle, colors)
  }
}

#' Draw a Single Clock Plot
#' @noRd
drawSingleClockPlot <- function(clockData, plotTitle, colors) {
  ggplot2::ggplot(clockData, ggplot2::aes(x = .data$Metric, y = .data$Value, fill = .data$Method)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::ylim(0, 1.1) +
    ggplot2::ggtitle(plotTitle) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "dashed"),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = .data$Value + 0.03, label = sprintf("%.3f", .data$Value)),
      position = ggplot2::position_dodge(width = 0.9), vjust = 0, size = 3
    )
}

#' Draw Multiple Clock Plots
#' @noRd
drawMultipleClockPlots <- function(clockData, plotTitle, colors) {
  methods <- unique(clockData$Method)
  if (length(methods) == 0) stop("No methods available to plot.")
  plotList <- lapply(methods, function(method) {
    methodData <- clockData[clockData$Method == method, ]
    # ZMIANA: Dodajemy prefiks stringr::
    methodData$Metric <- stringr::str_wrap(methodData$Metric, width = 10)
    titleText <- paste(plotTitle, "\n", method)
    # ZMIANA: Dodajemy prefiksy ggplot2::
    ggplot2::ggplot(methodData, ggplot2::aes(x = .data$Metric, y = .data$Value)) +
      ggplot2::geom_bar(stat = "identity", fill = ggplot2::alpha(colors[which(methods == method)], 0.5)) +
      ggplot2::coord_polar(start = 0) +
      ggplot2::ylim(0, 1.3) +
      ggplot2::ggtitle(titleText) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 20)),
        plot.margin = ggplot2::unit(c(1, 1, 2, 1), "cm")
      ) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Metric, y = 1.25), size = 3.5) +
      ggplot2::geom_text(ggplot2::aes(y = .data$Value + 0.05, label = sprintf("%.4f", .data$Value)), size = 3)
  })
  patchwork::wrap_plots(plotList)
}

#' Extract Metric Value from Confusion Matrix Object
#' @noRd
extractMetricValue <- function(cmObj, metricName) {
  if (metricName %in% names(cmObj$overall)) {
    cmObj$overall[[metricName]]
  } else if (metricName %in% names(cmObj$byClass)) {
    cmObj$byClass[[metricName]]
  } else {
    stop(sprintf("Metric '%s' not found in confusion matrix object.", metricName))
  }
}
