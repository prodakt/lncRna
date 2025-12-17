#' Create Clock Plots for Confusion Matrix Metrics
#'
#' Generates circular bar plots (clock plots) to visualize performance metrics.
#' Supports both single and multiple plot layouts for comparing methods.
#'
#' @param cmList A named list where each element represents a confusion matrix,
#'   typically the output from `calculateCM`.
#' @param methods An optional character vector of method names to include.
#'   If `NULL`, interactive selection is triggered in an interactive session;
#'   otherwise, all methods are used.
#' @param metrics A character vector of metrics to display. Defaults to
#'   `c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall")`.
#' @param plotTitle The title for the plot(s).
#' @param colors An optional vector of colors for each method.
#' @param layout The layout of the plots: `"single"` (default) for one plot,
#'   or `"multiple"` for a grid of plots.
#' @param ... Additional arguments (not currently used).
#'
#' @return A `ggplot` object (for `layout = "single"`) or a `patchwork`
#'   object (for `layout = "multiple"`), which can be printed to display.
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar ylim ggtitle
#' @importFrom ggplot2 scale_fill_manual theme_minimal theme element_blank
#' @importFrom ggplot2 element_line element_text geom_text position_dodge
#' @importFrom ggplot2 margin unit
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom scales alpha
#' @importFrom stringr str_wrap
#'
#' @examples
#' # --- 1. Create a mock cmList object (as from calculateCM) ---
#' set.seed(123)
#' mockCmList <- list(
#'   `MethodA` = list(metrics = c(Accuracy=0.9, Sensitivity=0.8, Specificity=0.95,
#'                                Precision=0.85, Recall=0.8)),
#'   `MethodB` = list(metrics = c(Accuracy=0.8, Sensitivity=0.9, Specificity=0.7,
#'                                Precision=0.75, Recall=0.9))
#' )
#'
#' # --- 2. Run the plot function ---
#' # Example 1: Single plot comparing selected methods
#' plotClockMetrics(
#'   cmList = mockCmList,
#'   methods = c("MethodA", "MethodB"),
#'   plotTitle = "Comparison Clock Plot"
#' )
#'
#' # Example 2: Multiple plots, one for each method
#' plotClockMetrics(cmList = mockCmList, layout = "multiple")
#'
plotClockMetrics <- function(cmList, methods = NULL,
                             metrics = c("Accuracy", "Sensitivity", "Specificity",
                                         "Precision", "Recall"),
                             plotTitle = "Clock Plot of Metrics", colors = NULL,
                             layout = "single", ...) {
    
    layout <- match.arg(layout, choices = c("single", "multiple"))
    
    validatedInputs <- validateClockInputs(cmList, methods, metrics)
    cmList <- validatedInputs$cmList
    methods <- validatedInputs$methods
    metrics <- validatedInputs$metrics
    
    if (length(methods) == 0) return(invisible(NULL))
    
    clockData <- prepareClockData(cmList, methods, metrics)
    
    plotColors <- setColors(colors, length(methods))
    
    if (layout == "single") {
        return(drawSingleClockPlot(clockData, plotTitle, plotColors, ...))
    } else {
        return(drawMultipleClockPlots(clockData, plotTitle, plotColors, ...))
    }
}

# --- Internal Helper Functions for Clock Plot ---

#' @noRd
validateClockInputs <- function(cmList, methods, metrics) {
    if (!is.list(cmList) || is.null(names(cmList)) || length(names(cmList)) == 0) {
        stop("'cmList' must be a non-empty, named list.")
    }
    
    if (!is.list(cmList[[1]]) || !("metrics" %in% names(cmList[[1]]))) {
        stop("Each element in 'cmList' must be a list containing a '$metrics' vector.")
    }
    
    availableMethods <- names(cmList)
    
    selectedMethods <- handleToolSelection(
        availableTools = availableMethods,
        selectedTools = methods
    )
    
    allAvailableMetrics <- names(cmList[[1]]$metrics)
    invalidMetrics <- setdiff(metrics, allAvailableMetrics)
    if (length(invalidMetrics) > 0) {
        stop("Invalid metrics specified: ", paste(invalidMetrics, collapse = ", "),
             ". Available metrics are: ", paste(allAvailableMetrics, collapse = ", "))
    }
    
    list(cmList = cmList, methods = selectedMethods, metrics = metrics)
}

#' @noRd
prepareClockData <- function(cmList, methods, metrics) {
    clockDataList <- lapply(methods, function(method) {
        values <- vapply(metrics, function(metric) {
            extractMetricValue(cmList[[method]], metric)
        }, FUN.VALUE = numeric(1))
        data.frame(
            Method = factor(method, levels = methods),
            Metric = factor(metrics, levels = metrics),
            Value = as.numeric(values)
        )
    })
    do.call(rbind, clockDataList)
}

#' @noRd
drawSingleClockPlot <- function(clockData, plotTitle, colors, ...) {
    Metric <- Value <- Method <- NULL
    ggplot2::ggplot(clockData, ggplot2::aes(x = Metric, y = Value, fill = Method)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        ggplot2::coord_polar(start = 0) +
        ggplot2::ylim(0, 1.1) +
        ggplot2::ggtitle(plotTitle) +
        ggplot2::scale_fill_manual(values = colors) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "dashed"),
            axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10, color = "black"),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        ) +
        ggplot2::geom_text(
            ggplot2::aes(y = Value + 0.05, label = sprintf("%.3f", Value)),
            position = ggplot2::position_dodge(width = 0.9),
            size = 3
        )
}

#' @noRd
drawMultipleClockPlots <- function(clockData, plotTitle, colors, ...) {
    Metric <- Value <- Method <- NULL
    methods <- levels(clockData$Method)
    
    plotList <- lapply(seq_along(methods), function(i) {
        methodName <- methods[i]
        methodData <- clockData[clockData$Method == methodName, ]
        
        ggplot2::ggplot(methodData, ggplot2::aes(x = Metric, y = Value)) +
            ggplot2::geom_bar(stat = "identity", fill = scales::alpha(colors[i], 0.7)) +
            ggplot2::coord_polar(start = 0) +
            ggplot2::ylim(0, 1.2) +
            ggplot2::ggtitle(paste(plotTitle, "-", methodName)) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "dashed"),
                axis.text.y = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(hjust = 0.5, size = 10)
            ) +
            ggplot2::geom_text(
                ggplot2::aes(y = Value + 0.05, label = sprintf("%.3f", Value)),
                size = 3
            )
    })
    
    patchwork::wrap_plots(plotList)
}