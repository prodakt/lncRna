#' Create Radar Plots for Confusion Matrix Metrics
#'
#' Generates radar plots (spider plots) to visualize and compare performance
#' metrics from confusion matrices. Requires the 'fmsb' and 'scales' packages.
#'
#' @param cmList A named list of `confusionMatrix` objects. Names of the list
#'   elements correspond to the method names.
#' @param methods An optional character vector of method names to plot. If `NULL`
#'   (default) and the session is interactive, the user is prompted to select
#'   methods. In a non-interactive session, `NULL` results in all methods
#'   being plotted.
#' @param metrics A character vector of metrics to visualize. Defaults to
#'   c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall").
#' @param plotTitle An optional title for the radar plot.
#' @param colors An optional character vector of colors for each method. If
#'   `NULL`, a default palette is used.
#' @param layout Specifies the layout: "single" for a combined plot or
#'   "multiple" for one plot per method. Defaults to "single".
#' @param displayArea A logical value; if `TRUE`, displays the polygon area.
#' @param displayFill A logical value; if `TRUE`, fills radar polygons with a
#'   semi-transparent color.
#' @param saveData A logical value; if `TRUE`, saves the radar data to a CSV.
#' @param fileName A character string specifying the file name if `saveData` is
#'   `TRUE`.
#' @param ... Additional arguments passed to `fmsb::radarchart`.
#'
#' @return The function is called for its side effect of creating a plot and
#'   returns `invisible(NULL)`. If required packages are missing, it also
#'   returns `invisible(NULL)` after printing a message.
#'
#' @export
#' @examples
#' if (requireNamespace("caret", quietly = TRUE) &&
#'     requireNamespace("fmsb", quietly = TRUE) &&
#'     requireNamespace("scales", quietly = TRUE)) {
#'
#'   # --- 1. Create a sample list of confusionMatrix objects ---
#'   set.seed(123)
#'   ref <- factor(sample(c("A", "B"), 50, replace = TRUE))
#'   pred1 <- factor(sample(c("A", "B"), 50, replace = TRUE))
#'   pred2 <- factor(sample(c("A", "B"), 50, replace = TRUE))
#'   sampleCmList <- list(Method1 = caret::confusionMatrix(pred1, ref),
#'                        Method2 = caret::confusionMatrix(pred2, ref))
#'
#'   # --- 2. Run non-interactive examples ---
#'   # Redirect graphics to a null device to prevent plots opening during checks
#'   png(tempfile())
#'
#'   radarPlotCm(
#'     cmList = sampleCmList,
#'     methods = "Method1",
#'     plotTitle = "Performance of Method1"
#'   )
#'
#'   dev.off() # Close the null device
#'
#' } else {
#'   message("Skipping radar plot examples because 'fmsb' or 'scales' is not installed.")
#' }
radarPlotCm <- function(cmList, methods = NULL,
                        metrics = c("Accuracy", "Sensitivity", "Specificity", "Precision", "Recall"),
                        plotTitle = "Radar Plot of Metrics", colors = NULL, layout = "single",
                        displayArea = FALSE, displayFill = TRUE, saveData = FALSE, fileName = NULL, ...) {

  # --- Check for suggested packages ---
  if (!requireNamespace("fmsb", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
    message("Packages 'fmsb' and 'scales' are required for the radarPlotCm() function.")
    message("Please install them using: BiocManager::install(c('fmsb', 'scales'))")
    return(invisible(NULL))
  }

  layout <- match.arg(layout, choices = c("single", "multiple"))
  validateRadarPlotInputs(cmList, metrics)
  availableMethods <- names(cmList)

  if (is.null(methods)) {
    selectedMethods <- selectRadarMethodsInteractively(availableMethods)
  } else {
    stopifnot("Some 'methods' are not in 'cmList'" = all(methods %in% availableMethods))
    selectedMethods <- methods
  }
  if (is.null(selectedMethods) || length(selectedMethods) == 0) {
    warning("No valid methods selected. Aborting plot.")
    return(invisible(NULL))
  }

  radarData <- prepareRadarData(cmList, selectedMethods, metrics)

  if (saveData) {
    stopifnot("'fileName' must be a character string when saveData is TRUE." =
                (!is.null(fileName) && is.character(fileName)))
    utils::write.csv2(radarData, file = fileName, row.names = TRUE)
    message("Radar data saved to ", fileName)
  }

  plotColors <- setRadarColors(colors, length(selectedMethods))
  dispatchRadarPlot(radarData, plotTitle, plotColors, layout, displayArea, displayFill, ...)
  return(invisible(NULL))
}

# --- Helper Functions ---

#' @noRd
selectRadarMethodsInteractively <- function(availableMethods) {
  if (!interactive()) {
    return(availableMethods)
  }
  message("Available methods:")
  choices <- c("All methods", availableMethods)
  choicesText <- paste0(seq_along(choices), ". ", choices)
  message(paste(choicesText, collapse = "\n"))
  prompt <- "Select method numbers (comma-separated): "
  userInput <- readline(prompt = prompt)
  if (userInput == "") {
    warning("No selection made.")
    return(NULL)
  }
  selectedIndices <- suppressWarnings(as.numeric(strsplit(userInput, ",")[[1]]))
  if (any(is.na(selectedIndices))) {
    warning("Invalid input. Please enter numbers.")
    return(NULL)
  }
  if (1 %in% selectedIndices) {
    return(availableMethods)
  }
  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(choices)]
  if (length(validIndices) == 0) {
    warning("Selected numbers do not correspond to available methods.")
    return(NULL)
  }
  return(availableMethods[validIndices - 1])
}

#' @noRd
validateRadarPlotInputs <- function(cmList, metrics) {
  stopifnot(
    "'cmList' must be a named list of confusionMatrix objects." =
      (is.list(cmList) && !is.null(names(cmList)) && length(names(cmList)) > 0),
    "All elements in 'cmList' must be of class 'confusionMatrix'." =
      all(sapply(cmList, inherits, "confusionMatrix"))
  )
  validMetrics <- unique(c("Accuracy", "Sensitivity", "Specificity", "Pos Pred Value",
                           "Neg Pred Value", "Precision", "Recall", "F1", "Prevalence",
                           "Detection Rate", "Detection Prevalence", "Balanced Accuracy",
                           names(cmList[[1]]$overall), names(cmList[[1]]$byClass)))
  stopifnot(
    "All specified 'metrics' must be valid." = all(metrics %in% validMetrics)
  )
}

#' @noRd
prepareRadarData <- function(cmList, methods, metrics) {
  radarDf <- data.frame(row.names = c("Max", "Min", methods))
  for (metric in metrics) {
    values <- sapply(methods, function(m) extractMetricValue(cmList[[m]], metric))
    if (any(values < 0 | values > 1, na.rm = TRUE)) {
      warning(sprintf("Metric '%s' has values outside [0,1]; scaling to fit.", metric))
      values <- scales::rescale(values, to = c(0, 1))
    }
    radarDf[[metric]] <- c(1, 0, as.numeric(values))
  }
  radarDf
}

#' @noRd
setRadarColors <- function(colors, nMethods) {
  if (is.null(colors)) {
    return(scales::hue_pal()(nMethods))
  }
  if (length(colors) < nMethods) {
    warning("Fewer colors than methods; recycling colors.")
    return(rep_len(colors, nMethods))
  }
  return(colors[seq_len(nMethods)])
}

#' @noRd
dispatchRadarPlot <- function(radarData, plotTitle, colors, layout, displayArea, displayFill, ...) {
  if (layout == "single") {
    drawSingleRadarPlot(radarData, plotTitle, colors, displayArea, displayFill, ...)
  } else {
    drawMultipleRadarPlots(radarData, plotTitle, colors, displayArea, displayFill, ...)
  }
}

#' @noRd
drawSingleRadarPlot <- function(radarData, plotTitle, colors, displayArea, displayFill, ...) {
  plotFill <- if (displayFill) scales::alpha(colors, 0.5) else NA
  fmsb::radarchart(radarData, pfcol = plotFill, pcol = colors, plwd = 2, cglcol = "grey",
                   title = plotTitle, ...)

  legendText <- rownames(radarData)[-c(1, 2)]
  if (displayArea) {
    areas <- sapply(legendText, function(method) {
      areaVal <- calculateRadarPolygonArea(radarData, method, colnames(radarData))
      paste0(method, " (Area=", sprintf("%.4f", areaVal), ")")
    })
    legendText <- areas
  }
  graphics::legend("bottom", legend = legendText, bty = "n", pch = 20, col = colors,
                   text.col = "black", cex = 0.8, pt.cex = 1.5, xpd = TRUE)
}

#' @noRd
drawMultipleRadarPlots <- function(radarData, plotTitle, colors, displayArea, displayFill, ...) {
  numMethods <- nrow(radarData) - 2
  stopifnot("No methods available to plot." = (numMethods > 0))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

  layoutDims <- calculateOptimalLayout(numMethods)
  graphics::par(mar = c(3, 3, 3, 1), mfrow = layoutDims)

  methodsToPlot <- rownames(radarData)[-c(1, 2)]
  for (i in seq_along(methodsToPlot)) {
    methodName <- methodsToPlot[i]
    methodData <- radarData[c(1, 2, which(rownames(radarData) == methodName)), ]
    plotFill <- if (displayFill) scales::alpha(colors[i], 0.5) else NA

    fmsb::radarchart(methodData, pfcol = plotFill, pcol = colors[i], plwd = 2,
                     title = paste(plotTitle, "-", methodName), ...)

    if (displayArea) {
      polygonArea <- calculateRadarPolygonArea(radarData, methodName, colnames(radarData))
      graphics::mtext(paste("Area =", sprintf("%.4f", polygonArea)), side = 1, line = 1,
                      cex = 0.8, col = colors[i])
    }
  }
}

#' @noRd
calculateRadarPolygonArea <- function(radarData, methodName, metrics) {
  methodValues <- as.numeric(radarData[rownames(radarData) == methodName, metrics])
  nMetrics <- length(methodValues)
  angles <- seq(0, 2 * pi, length.out = nMetrics + 1)[seq_len(nMetrics)]
  xCoords <- methodValues * cos(angles)
  yCoords <- methodValues * sin(angles)
  area <- 0.5 * abs(sum(xCoords * c(yCoords[-1], yCoords[1])) -
                      sum(yCoords * c(xCoords[-1], xCoords[1])))
  area
}

#' @noRd
calculateOptimalLayout <- function(nPlots) {
  if (nPlots <= 2) return(c(1, nPlots))
  if (nPlots <= 4) return(c(2, 2))
  if (nPlots <= 6) return(c(2, 3))
  rows <- ceiling(sqrt(nPlots))
  cols <- ceiling(nPlots / rows)
  c(rows, cols)
}

#' @noRd
extractMetricValue <- function(cmObj, metricName) {
  if (metricName %in% names(cmObj$overall)) {
    cmObj$overall[[metricName]]
  } else if (metricName %in% names(cmObj$byClass)) {
    cmObj$byClass[[metricName]]
  } else {
    stop(sprintf("Metric '%s' not found in confusion matrix.", metricName))
  }
}
