#' Plot Sankey Diagram for Selected lncRNAs
#'
#' Creates a Sankey diagram visualizing interactions for specific lncRNAs from
#' a gProfiler output table. Allows highlighting or exclusive plotting of
#' selected lncRNAs.
#'
#' @param data A data.frame containing gProfiler interaction results. Must
#'   include columns: `term_name`, `intersection`, `lncRNAId`, and `type`.
#' @param selectedLnc A character vector of lncRNA IDs to highlight or plot.
#'   If `NULL` (default), all interactions are plotted.
#' @param showLabel A logical value; if `TRUE`, shows node labels. Defaults to
#'   `FALSE`.
#' @param nodeColor A character string specifying a single color for all nodes.
#'   If `NULL` (default), a random palette is used.
#' @param title A character string for the diagram title. Defaults to `NULL`.
#' @param colorSelected A logical value. If `TRUE`, all interactions are plotted,
#'   but only those related to `selectedLnc` are colored (others are gray). If
#'   `FALSE` (default), only interactions involving `selectedLnc` are plotted.
#'   Requires `selectedLnc` to be specified.
#'
#' @return A `plotly` object representing the Sankey diagram.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom Polychrome createPalette
#' @export
#' @examples
#' # --- 1. Create a sample data frame ---
#' sampleData <- data.frame(
#'   term_name = c("GO:1", "GO:2", "GO:3", "GO:2"),
#'   intersection = c("GeneA", "GeneB", "GeneC", "GeneC"),
#'   lncRNAId = c("Lnc1", "Lnc1", "Lnc2", "Lnc3"),
#'   type = c("cis", "trans", "cis", "trans"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Generate Sankey plots ---
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("Polychrome", quietly = TRUE)) {
#'
#'   # Example 2a: Plot only interactions for Lnc1 and Lnc2
#'   fig1 <- plotByLnc(
#'     data = sampleData,
#'     selectedLnc = c("Lnc1", "Lnc2"),
#'     title = "Interactions for Lnc1 & Lnc2"
#'   )
#'
#'   # Example 2b: Plot all interactions, but highlight those for Lnc1
#'   fig2 <- plotByLnc(
#'     data = sampleData,
#'     selectedLnc = "Lnc1",
#'     colorSelected = TRUE,
#'     title = "All Interactions (Lnc1 Highlighted)"
#'   )
#'
#'   # In a real session, you would just type `fig1` or `fig2` to display.
#' }
plotByLnc <- function(data, selectedLnc = NULL, showLabel = FALSE,
                      nodeColor = NULL, title = NULL, colorSelected = FALSE) {

  # --- 1. Validate inputs and prepare data ---
  preparedData <- validateLncPlotInputs(data, selectedLnc, colorSelected)
  if (is.null(preparedData)) {
    return(invisible(NULL))
  }

  # --- 2. Prepare Sankey data ---
  sankeyData <- prepareLncPlotData(preparedData$data)

  # --- 3. Set colors for nodes and links ---
  sankeyData <- setLncPlotSankeyColors(
    sankeyData, nodeColor, colorSelected, preparedData$pathway
  )

  # --- 4. Create the Sankey plot ---
  createLncPlotSankey(sankeyData, showLabel, title)
}

# --- Helper Functions ---

#' Validate and Prepare Data for lncRNA Plotting
#' @noRd
validateLncPlotInputs <- function(data, selectedLnc, colorSelected) {
  requiredCols <- c("intersection", "term_name", "type", "lncRNAId")
  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' must contain required columns" = all(requiredCols %in% colnames(data))
  )

  if (colorSelected && is.null(selectedLnc)) {
    stop("To use 'colorSelected = TRUE', please provide 'selectedLnc'.")
  }

  pathwayData <- NULL

  if (!is.null(selectedLnc)) {
    if (!any(selectedLnc %in% data$lncRNAId)) {
      warning("None of the selected lncRNAs were found in the data.")
      return(NULL)
    }

    if (colorSelected) {
      # Keep all data, but define the pathway for coloring
      pathwayData <- data[data$lncRNAId %in% selectedLnc, ]
    } else {
      # Filter data to only include selected lncRNAs
      data <- data[data$lncRNAId %in% selectedLnc, ]
    }
  }

  return(list(data = data, pathway = pathwayData))
}

#' Prepare Data for lncRNA Sankey Plot
#' @noRd
prepareLncPlotData <- function(data) {
  goProt <- data.frame(source = data$intersection, target = data$term_name)
  lncProt <- data.frame(source = data$lncRNAId, target = data$intersection)

  link <- rbind(goProt, lncProt)
  link$value <- 1

  nodeNames <- unique(c(link$source, link$target))
  node <- data.frame(name = nodeNames)

  link$IDsource <- match(link$source, nodeNames) - 1
  link$IDtarget <- match(link$target, nodeNames) - 1

  list(link = link, node = node)
}

#' Set Colors for lncRNA Sankey Plot
#' @noRd
setLncPlotSankeyColors <- function(sankeyData, nodeColor, colorSelected, pathway) {
  node <- sankeyData$node
  link <- sankeyData$link

  if (is.null(nodeColor)) {
    node$color <- Polychrome::createPalette(nrow(node), c("#ff0000", "#00ff00", "#0000ff"))
  } else {
    node$color <- nodeColor
  }

  if (colorSelected && !is.null(pathway)) {
    selectedNodes <- unique(c(pathway$intersection, pathway$term_name, pathway$lncRNAId))
    node$color[!node$name %in% selectedNodes] <- "gray"
  }

  link$color <- node$color[link$IDtarget + 1]

  if (colorSelected && !is.null(pathway)) {
    # This logic might need refinement based on desired coloring effect.
    # Here, we gray out links not involving the selected pathway's components.
    selectedSources <- unique(c(pathway$lncRNAId, pathway$intersection))
    selectedTargets <- unique(c(pathway$intersection, pathway$term_name))
    isLinkInPathway <- (link$source %in% selectedSources) & (link$target %in% selectedTargets)
    link$color[!isLinkInPathway] <- "gray"
  }

  list(node = node, link = link)
}

#' Create lncRNA Sankey Plot with Plotly
#' @noRd
createLncPlotSankey <- function(sankeyData, showLabel, title) {
  plotNode <- list(
    label = if (showLabel) sankeyData$node$name else NULL,
    color = sankeyData$node$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5)
  )
  plotLink <- list(
    source = sankeyData$link$IDsource,
    target = sankeyData$link$IDtarget,
    value = sankeyData$link$value,
    color = sankeyData$link$color
  )

  plotly::plot_ly(
    type = "sankey",
    orientation = "h",
    node = plotNode,
    link = plotLink
  ) |>
    plotly::layout(
      title = title,
      font = list(size = 10)
    )
}
