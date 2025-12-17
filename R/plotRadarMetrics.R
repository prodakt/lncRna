#' Create Radar Plots for Confusion Matrix Metrics
#'
#' Generates radar plots to visualize and compare performance metrics for
#' different methods or tool combinations. Supports a single plot for comparing
#' multiple methods or a grid of plots for individual evaluation.
#'
#' @param cmList A named list where each element represents a confusion matrix,
#'   typically the output from `calculateCM`. Each element must be a list
#'   containing at least a named numeric vector called `$metrics`.
#' @param methods An optional character vector of method names (from `names(cmList)`)
#'   to include in the plot. If `NULL`, behavior depends on the session:
#'   interactive mode prompts for selection, while non-interactive mode uses
#'   all available methods.
#' @param metrics A character vector of metrics to visualize on the radar plot axes.
#'   Defaults to `c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall")`.
#' @param plotTitle A character string for the plot title.
#' @param colors An optional character vector of colors for each method. If `NULL`,
#'   default colors are generated.
#' @param layout Specifies the plot layout: `"single"` (default) for one plot
#'   comparing all methods, or `"multiple"` for a grid of individual plots.
#' @param displayArea Logical. If `TRUE`, displays the calculated polygon area.
#' @param displayFill Logical. If `TRUE` (default), fills radar chart polygons
#'   with a semi-transparent color.
#' @param saveData Logical. If `TRUE`, saves the underlying data frame to a CSV file.
#' @param fileName The name of the file for saving data (required if `saveData` is `TRUE`).
#' @param ... Additional arguments passed to `fmsb::radarchart`.
#'
#' @return This function is called for its side effect of generating plots and
#'   does not return a value (`invisible(NULL)`).
#' @export
#' @importFrom fmsb radarchart
#' @importFrom scales alpha
#' @importFrom graphics legend mtext par
#' @importFrom utils write.csv2
#'
#' @examples
#' # --- 1. Create a mock cmList object (as from calculateCM) ---
#' set.seed(123)
#' mockCmList <- list(
#'   `MethodA` = list(metrics = c(Accuracy=0.9, Sensitivity=0.8, Specificity=0.95,
#'                                Precision=0.85, Recall=0.8)),
#'   `MethodB` = list(metrics = c(Accuracy=0.8, Sensitivity=0.9, Specificity=0.7,
#'                                Precision=0.75, Recall=0.9)),
#'   `MethodC` = list(metrics = c(Accuracy=0.85, Sensitivity=0.85, Specificity=0.85,
#'                                Precision=0.85, Recall=0.85))
#' )
#'
#' # --- 2. Run the plot function ---
#' # To prevent plots from showing up during automated checks, we wrap in a device
#' temp_png <- tempfile(fileext = ".png")
#' png(temp_png)
#'
#' # Example 1: Single plot comparing selected methods
#' plotRadarMetrics(
#'   cmList = mockCmList,
#'   methods = c("MethodA", "MethodC"),
#'   plotTitle = "Comparison Plot"
#' )
#'
#' # Example 2: Multiple plots, one for each method
#' plotRadarMetrics(cmList = mockCmList, layout = "multiple")
#'
#' dev.off()
#' unlink(temp_png)
#'
plotRadarMetrics <- function(cmList, methods = NULL,
                             metrics = c("Accuracy", "Sensitivity", "Specificity",
                                         "Precision", "Recall"),
                             plotTitle = "Radar Plot of Metrics", colors = NULL,
                             layout = "single", displayArea = FALSE,
                             displayFill = TRUE, saveData = FALSE,
                             fileName = NULL, ...) {
    
    layout <- match.arg(layout, choices = c("single", "multiple"))
    
    validatedInputs <- validateRadarInputs(cmList, methods, metrics)
    methods <- validatedInputs$methods
    metrics <- validatedInputs$metrics
    
    if (length(methods) == 0) return(invisible(NULL))
    
    radarData <- prepareRadarData(cmList, methods, metrics)
    
    if (saveData) {
        if (is.null(fileName) || !is.character(fileName)) {
            stop("A valid file name is required when saveData is TRUE.")
        }
        utils::write.csv2(radarData, file = paste0(fileName, ".csv"))
        message("Radar data frame saved to ", paste0(fileName, ".csv"))
    }
    
    plotColors <- setColors(colors, length(methods))
    
    drawRadarPlot(radarData, plotTitle, plotColors, layout, displayArea,
                  displayFill, ...)
    
    return(invisible(NULL))
}


# --- Internal Helper Functions Specific to Radar Plot ---

#' @noRd
validateRadarInputs <- function(cmList, methods, metrics) {
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
prepareRadarData <- function(cmList, methods, metrics) {
    radar_df <- data.frame(row.names = c("Max", "Min", methods))
    
    for (metric in metrics) {
        values <- vapply(methods, function(m) {
            extractMetricValue(cmList[[m]], metric)
        }, FUN.VALUE = numeric(1))
        
        if (any(values < 0 | values > 1, na.rm = TRUE)) {
            warning("Metric '", metric, "' has values outside [0,1]; scaling to fit.",
                    call. = FALSE)
            values <- scales::rescale(values, to = c(0, 1), from = range(values, na.rm = TRUE))
        }
        
        radar_df[[metric]] <- c(1, 0, values)
    }
    return(radar_df)
}

#' @noRd
drawRadarPlot <- function(radarData, plotTitle, colors, layout,
                          displayArea, displayFill, ...) {
    if (layout == "single") {
        drawSingleRadar(radarData, plotTitle, colors, displayArea, displayFill, ...)
    } else {
        drawMultipleRadars(radarData, plotTitle, colors, displayArea, displayFill, ...)
    }
}

#' @noRd
drawSingleRadar <- function(radarData, plotTitle, colors, displayArea, displayFill, ...) {
    plotFill <- if (displayFill) scales::alpha(colors, 0.5) else NA
    
    fmsb::radarchart(
        radarData, pfcol = plotFill, pcol = colors, plwd = 2, plty = 1,
        cglcol = "grey", cglty = 1, axislabcol = "black",
        caxislabels = seq(0, 1, 0.2), cglwd = 0.8, axistype = 1,
        seg = 5, title = plotTitle, ...
    )
    
    legendText <- rownames(radarData)[-c(1, 2)]
    if (displayArea) {
        legendText <- vapply(legendText, function(method) {
            areaVal <- calculatePolygonArea(radarData, method)
            paste0(method, " (Area=", sprintf("%.4f", areaVal), ")")
        }, FUN.VALUE = character(1))
    }
    
    graphics::legend("bottom", legend = legendText, bty = "n", pch = 20, col = colors,
                     text.col = "black", cex = 0.8, pt.cex = 1.2, xpd = TRUE)
}

#' @noRd
drawMultipleRadars <- function(radarData, plotTitle, colors, displayArea, displayFill, ...) {
    numMethods <- nrow(radarData) - 2
    if (numMethods <= 0) return()
    
    methodsToPlot <- rownames(radarData)[-c(1, 2)]
    
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    
    num_plots_per_page <- prod(calculateOptimalLayout(min(numMethods, 9)))
    
    graphics::par(mfrow = calculateOptimalLayout(min(numMethods, num_plots_per_page)))
    
    for (i in seq_along(methodsToPlot)) {
        if (i > 1 && (i - 1) %% num_plots_per_page == 0) {
            if (interactive()) {
                readline(prompt = "Press [Enter] to see the next page of plots...")
            }
            remaining_plots <- numMethods - (i - 1)
            graphics::par(mfrow = calculateOptimalLayout(min(remaining_plots, num_plots_per_page)))
        }
        
        methodName <- methodsToPlot[i]
        methodData <- radarData[c(1, 2, i + 2), ]
        plotFill <- if (displayFill) scales::alpha(colors[i], 0.5) else NA
        
        fmsb::radarchart(
            methodData, pfcol = plotFill, pcol = colors[i], plwd = 2, plty = 1,
            cglcol = "grey", cglty = 1, axislabcol = "grey",
            caxislabels = seq(0, 1, 0.25), cglwd = 0.8, axistype = 1,
            seg = 4, title = paste(plotTitle, "\n-", methodName), ...
        )
        
        if (displayArea) {
            areaVal <- calculatePolygonArea(radarData, methodName)
            graphics::mtext(paste("Area =", sprintf("%.4f", areaVal)), side = 1,
                            line = -1, cex = 0.8, col = colors[i])
        }
    }
}

#' @noRd
calculatePolygonArea <- function(radarData, methodName) {
    methodValues <- as.numeric(radarData[methodName, ])
    nMetrics <- length(methodValues)
    if (nMetrics < 3) return(NA_real_)
    
    angles <- seq(0, 2 * pi, length.out = nMetrics + 1)[- (nMetrics + 1)]
    xCoords <- methodValues * cos(angles)
    yCoords <- methodValues * sin(angles)
    
    area <- 0.5 * abs(sum(xCoords * c(yCoords[-1], yCoords[1])) -
                          sum(yCoords * c(xCoords[-1], xCoords[1])))
    return(area)
}