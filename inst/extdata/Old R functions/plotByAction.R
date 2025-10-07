#' Plot Sankey Diagram by Interaction Type
#'
#' Creates a Sankey diagram visualizing interactions from a gProfiler output
#' table, filtered by interaction type (cis/trans).
#'
#' @param data A data.frame containing gProfiler interaction results. Must
#'   include columns: `term_name`, `intersection`, `lncRNAId`, and `type`.
#' @param cis A logical value; if `TRUE` (default), includes cis interactions.
#' @param trans A logical value; if `TRUE` (default), includes trans interactions.
#' @param showLabel A logical value; if `TRUE`, shows node labels. Defaults to
#'   `FALSE`.
#' @param nodeColor A character string specifying a single color for all nodes
#'   and links. If `NULL` (default), a random palette is used.
#' @param title A character string for the diagram title. Defaults to `NULL`.
#'
#' @return A `plotly` object representing the Sankey diagram.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom Polychrome createPalette
#' @export
#' @examples
#' # --- 1. Create a sample data frame ---
#' sampleData <- data.frame(
#'   term_name = c("GO:123", "GO:456", "GO:123"),
#'   intersection = c("GeneA", "GeneB", "GeneB"),
#'   lncRNAId = c("Lnc1", "Lnc1", "Lnc2"),
#'   type = c("cis", "trans", "cis"),
#'   source = c("GO:BP", "GO:BP", "GO:BP"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Generate a Sankey plot for cis interactions only ---
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("Polychrome", quietly = TRUE)) {
#'   fig <- plotByAction(
#'     data = sampleData,
#'     cis = TRUE,
#'     trans = FALSE,
#'     title = "Cis-Regulatory Interactions"
#'   )
#'   # In a real session, you would just type `fig` to display the plot.
#'   # print(fig) # Printing may be needed in some contexts.
#' }
#'
plotByAction <- function(data, cis = TRUE, trans = TRUE, showLabel = FALSE,
                         nodeColor = NULL, title = NULL) {
  # --- 1. Validate inputs and filter data ---
  data <- validateActionPlotInputs(data, cis, trans)
  if (is.null(data)) {
    # Stop if validation returns NULL (e.g., no data left)
    return(invisible(NULL))
  }

  # --- 2. Prepare Sankey data ---
  sankeyData <- prepareActionPlotData(data)

  # --- 3. Set colors for nodes and links ---
  sankeyData <- setActionPlotSankeyColors(sankeyData, nodeColor)

  # --- 4. Create the Sankey plot ---
  createActionPlotSankey(sankeyData, showLabel, title)
}

# --- Helper Functions ---

#' Validate Inputs for Plotting by Interaction Type
#' @noRd
validateActionPlotInputs <- function(data, cis, trans) {
  requiredCols <- c("intersection", "term_name", "type", "lncRNAId")
  if (!all(requiredCols %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(requiredCols, collapse = ", "))
  }
  if (!cis && !trans) {
    stop("At least one of 'cis' or 'trans' must be TRUE.")
  }

  typesToKeep <- c()
  if (cis) typesToKeep <- c(typesToKeep, "cis")
  if (trans) typesToKeep <- c(typesToKeep, "trans")

  data <- data[data$type %in% typesToKeep, ]

  if (nrow(data) == 0) {
    warning("No data remains after filtering for specified interaction types.")
    return(NULL)
  }
  return(data)
}

#' Prepare Data for Sankey Plot
#' @noRd
prepareActionPlotData <- function(data) {
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

#' Set Colors for Sankey Plot
#' @noRd
setActionPlotSankeyColors <- function(sankeyData, nodeColor) {
  if (is.null(nodeColor)) {
    sankeyData$node$color <- Polychrome::createPalette(
      nrow(sankeyData$node), c("#ff0000", "#00ff00", "#0000ff")
    )
  } else {
    sankeyData$node$color <- nodeColor
  }
  sankeyData$link$color <- sankeyData$node$color[sankeyData$link$IDtarget + 1]
  return(sankeyData)
}

#' Create Sankey Plot with Plotly
#' @noRd
createActionPlotSankey <- function(sankeyData, showLabel, title) {
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
