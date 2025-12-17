#' Plot Sankey Diagram for gProfiler Term Connections
#'
#' Creates a Sankey diagram visualizing connections to selected gProfiler terms.
#' Requires the 'plotly' and 'Polychrome' packages.
#'
#' @param data A data.frame from interaction processing functions. Must include
#'   `term_name`, `intersection`, `lncRNAId`, and `type` columns.
#' @param selectedTerms A character vector of term names to highlight or plot.
#'   If `NULL` (default) and the session is interactive, the user is prompted
#'   to select terms. In a non-interactive session, `NULL` results in all
#'   terms being plotted.
#' @param showLabel A logical value; if `TRUE`, node labels are displayed.
#' @param nodeColor A single color for all nodes. If `NULL` (default), a
#'   random palette is used.
#' @param title A character string for the diagram title.
#' @param colorSelected A logical value. If `TRUE`, all interactions are plotted,
#'   but only those related to `selectedTerms` are colored. If `FALSE`
#'   (default), only interactions involving `selectedTerms` are plotted.
#'
#' @return A `plotly` object representing the Sankey diagram, or
#'   `invisible(NULL)` if required packages are not installed.
#'
#' @export
#' @examples
#' # --- 1. Create a sample data frame for all examples ---
#' sampleData <- data.frame(
#'   term_name = c("GO:1", "GO:2", "GO:3", "GO:2"),
#'   intersection = c("GeneA", "GeneB", "GeneC", "GeneA"),
#'   lncRNAId = c("Lnc1", "Lnc1", "Lnc2", "Lnc2"),
#'   type = c("cis", "trans", "cis", "trans"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Demonstrate plotting functions (will run only if packages are available) ---
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("Polychrome", quietly = TRUE)) {
#'
#'   fig_terms <- plotByTerms(data = sampleData, selectedTerms = "GO:2")
#'   # To display the plot in an interactive session, simply type `fig_terms`
#'
#'   fig_target <- plotByTarget(data = sampleData, selectedTarget = "GeneA")
#'
#' } else {
#'   message("Skipping Sankey plot examples because 'plotly' or 'Polychrome' is not installed.")
#' }
plotByTerms <- function(data, selectedTerms = NULL, showLabel = FALSE,
                        nodeColor = NULL, title = NULL, colorSelected = FALSE) {
  if (!checkSankeyDependencies()) return(invisible(NULL))

  by_col <- "term_name"
  validateSankeyInputs(data, selectedTerms, colorSelected, by = by_col)

  if (is.null(selectedTerms)) {
    selectedTerms <- selectSankeyItemsInteractively(data, by = by_col)
    if (is.null(selectedTerms)) return(invisible(NULL))
  }

  preparedData <- prepareSankeyHighlightData(data, selectedTerms, colorSelected, by = by_col)
  if (is.null(preparedData)) return(invisible(NULL))

  createSankeyPlot(
    data = preparedData$data, pathwayData = preparedData$pathway,
    colorSelected = colorSelected, showLabel = showLabel,
    nodeColor = nodeColor, title = title, by = by_col
  )
}

#' @rdname plotByTerms
#' @param selectedTarget A character vector of target gene IDs to highlight or plot.
#' @export
plotByTarget <- function(data, selectedTarget = NULL, showLabel = FALSE,
                         nodeColor = NULL, title = NULL, colorSelected = FALSE) {
  if (!checkSankeyDependencies()) return(invisible(NULL))

  by_col <- "intersection"
  validateSankeyInputs(data, selectedTarget, colorSelected, by = by_col)

  if (is.null(selectedTarget)) {
    selectedTarget <- selectSankeyItemsInteractively(data, by = by_col)
    if (is.null(selectedTarget)) return(invisible(NULL))
  }

  preparedData <- prepareSankeyHighlightData(data, selectedTarget, colorSelected, by = by_col)
  if (is.null(preparedData)) return(invisible(NULL))

  createSankeyPlot(
    data = preparedData$data, pathwayData = preparedData$pathway,
    colorSelected = colorSelected, showLabel = showLabel,
    nodeColor = nodeColor, title = title, by = by_col
  )
}

#' @rdname plotByTerms
#' @param selectedLnc A character vector of lncRNA IDs to highlight or plot.
#' @export
plotByLnc <- function(data, selectedLnc = NULL, showLabel = FALSE,
                      nodeColor = NULL, title = NULL, colorSelected = FALSE) {
  if (!checkSankeyDependencies()) return(invisible(NULL))

  by_col <- "lncRNAId"
  validateSankeyInputs(data, selectedLnc, colorSelected, by = by_col)

  if (is.null(selectedLnc)) {
    selectedLnc <- selectSankeyItemsInteractively(data, by = by_col)
    if (is.null(selectedLnc)) return(invisible(NULL))
  }

  preparedData <- prepareSankeyHighlightData(data, selectedLnc, colorSelected, by = by_col)
  if (is.null(preparedData)) return(invisible(NULL))

  createSankeyPlot(
    data = preparedData$data, pathwayData = preparedData$pathway,
    colorSelected = colorSelected, showLabel = showLabel,
    nodeColor = nodeColor, title = title, by = by_col
  )
}

#' @rdname plotByTerms
#' @param interactionTypes A character vector specifying the interaction types to plot.
#' @export
plotByType <- function(data, interactionTypes = NULL, showLabel = FALSE,
                       nodeColor = NULL, title = NULL) {
  if (!checkSankeyDependencies()) return(invisible(NULL))

  by_col <- "type"
  validateSankeyInputs(data, interactionTypes, colorSelected = FALSE, by = by_col)

  if (is.null(interactionTypes)) {
    interactionTypes <- selectSankeyItemsInteractively(data, by = by_col)
    if (is.null(interactionTypes)) return(invisible(NULL))
  }

  preparedData <- prepareSankeyHighlightData(data, interactionTypes, colorSelected = FALSE, by = by_col)
  if (is.null(preparedData)) return(invisible(NULL))

  createSankeyPlot(
    data = preparedData$data, pathwayData = NULL,
    colorSelected = FALSE, showLabel = showLabel,
    nodeColor = nodeColor, title = title, by = by_col
  )
}


# --- Generic Helper Functions (Internal, for all Sankey plots) ---

#' Check for Sankey Plot Dependencies
#' @return A logical value: TRUE if all dependencies are met, FALSE otherwise.
#' @noRd
checkSankeyDependencies <- function() {
  if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("Polychrome", quietly = TRUE)) {
    message("Packages 'plotly' and 'Polychrome' are required for this plotting function.")
    message("Please install them using: BiocManager::install(c('plotly', 'Polychrome'))")
    return(FALSE)
  }
  return(TRUE)
}


#' Core Sankey Plotting Function (Internal)
#' @noRd
createSankeyPlot <- function(data, pathwayData, colorSelected, showLabel,
                             nodeColor, title, by) {
  goProt <- data.frame(source = data$intersection, target = data$term_name)
  lncProt <- data.frame(source = data$lncRNAId, target = data$intersection)
  link <- rbind(goProt, lncProt)
  link$value <- 1
  nodeNames <- unique(c(link$source, link$target))
  node <- data.frame(name = nodeNames)
  link$IDsource <- match(link$source, nodeNames) - 1
  link$IDtarget <- match(link$target, nodeNames) - 1

  if (is.null(nodeColor)) {
    nNodes <- nrow(node)
    palette <- Polychrome::createPalette(max(nNodes, 3), c("#ff0000", "#00ff00", "#0000ff"))
    node$color <- palette[seq_len(nNodes)]
  } else {
    node$color <- nodeColor
  }
  if (colorSelected && !is.null(pathwayData)) {
    pathwayNodes <- unique(c(pathwayData$intersection, pathwayData$term_name, pathwayData$lncRNAId))
    node$color[!node$name %in% pathwayNodes] <- "gray"
    link$color <- node$color[link$IDtarget + 1]
    highlightColumn <- if (by == "term_name") "term_name" else "intersection"
    link$color[!link$target %in% pathwayData[[highlightColumn]]] <- "gray"
  } else {
    link$color <- node$color[link$IDtarget + 1]
  }

  plotNode <- list(
    label = if (showLabel) node$name else NULL,
    color = node$color, pad = 15, thickness = 20,
    line = list(color = "black", width = 0.5)
  )
  plotLink <- list(
    source = link$IDsource, target = link$IDtarget,
    value = link$value, color = link$color
  )
  plotly::plot_ly(type = "sankey", orientation = "h", node = plotNode, link = plotLink) |>
    plotly::layout(title = title, font = list(size = 10))
}

#' Validate Inputs for Sankey Plot Functions (Internal)
#' @noRd
validateSankeyInputs <- function(data, selectedItems, colorSelected, by) {
  requiredCols <- c("intersection", "term_name", "type", "lncRNAId")
  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' must contain required columns" = all(requiredCols %in% colnames(data))
  )
  if (colorSelected && is.null(selectedItems) && by != "type") {
    stop("To use 'colorSelected = TRUE', please provide selected items.")
  }
  if (!is.null(selectedItems) && !any(selectedItems %in% data[[by]])) {
    warning("None of the selected items were found in the data.")
  }
}

#' Interactively Select Items or Default to All (Internal)
#' @noRd
selectSankeyItemsInteractively <- function(data, by) {
  availableItems <- unique(data[[by]])
  if (!interactive()) {
    return(availableItems)
  }
  itemLabel <- switch(by,
                      "term_name" = "terms",
                      "intersection" = "targets",
                      "lncRNAId" = "lncRNAs",
                      "type" = "types",
                      "items"
  )
  message("Select ", itemLabel, " to plot:")
  choices <- c(paste("All", itemLabel), availableItems)
  choicesText <- paste0(seq_along(choices), ". ", choices)
  message(paste(choicesText, collapse = "\n"))
  prompt <- "Enter numbers (comma-separated): "
  userInput <- readline(prompt = prompt)
  if (userInput == "") {
    warning("No selection made.")
    return(NULL)
  }
  selectedIndices <- suppressWarnings(as.integer(strsplit(userInput, ",")[[1]]))
  if (any(is.na(selectedIndices))) {
    warning("Invalid input. Please enter numbers.")
    return(NULL)
  }
  if (1 %in% selectedIndices) {
    return(availableItems)
  }
  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(choices)]
  if (length(validIndices) == 0) {
    warning("Selected numbers do not correspond to available items.")
    return(NULL)
  }
  return(availableItems[validIndices - 1])
}

#' Prepare Data for Filtering or Highlighting (Internal)
#' @noRd
prepareSankeyHighlightData <- function(data, selectedItems, colorSelected, by) {
  validSelectedItems <- intersect(selectedItems, data[[by]])
  if (length(validSelectedItems) == 0) {
    warning("None of the selected items are in the data.")
    return(NULL)
  }
  if (colorSelected) {
    pathway <- data[data[[by]] %in% validSelectedItems, ]
    return(list(data = data, pathway = pathway))
  } else {
    filteredData <- data[data[[by]] %in% validSelectedItems, ]
    return(list(data = filteredData, pathway = NULL))
  }
}
