#' Plot Sankey Diagram by Specified Interaction Types
#'
#' Creates a Sankey diagram visualizing interactions from a gProfiler output
#' table, filtered by user-specified interaction types.
#'
#' @param data A data.frame containing gProfiler interaction results. Must
#'   include columns: `term_name`, `intersection`, `lncRNAId`, and `type`.
#' @param interactionTypes A character vector specifying the interaction types
#'   to include (e.g., "cis", "trans"). If `NULL` (default) and the session is
#'   interactive, the user will be prompted to select types. In a
#'   non-interactive session, `NULL` results in all types being plotted.
#' @param showLabel A logical value; if `TRUE`, shows node labels. Defaults to
#'   `FALSE`.
#' @param nodeColor A single color for all nodes and links. If `NULL`
#'   (default), a random palette is used.
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
#'   term_name = c("GO:1", "GO:2", "GO:3"),
#'   intersection = c("GeneA", "GeneB", "GeneA"),
#'   lncRNAId = c("Lnc1", "Lnc1", "Lnc2"),
#'   type = c("cis", "trans", "cis"),
#'   stringsAsFactors = FALSE
#' )
#'
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("Polychrome", quietly = TRUE)) {
#'
#'   # --- Example 1: Non-interactive plot for a specific type ---
#'   fig1 <- plotByType(
#'     data = sampleData,
#'     interactionTypes = "cis",
#'     title = "Cis Interactions"
#'   )
#'
#'   # --- Example 2: Non-interactive plot for multiple types ---
#'   fig2 <- plotByType(
#'     data = sampleData,
#'     interactionTypes = c("cis", "trans"),
#'     showLabel = TRUE,
#'     title = "All Interaction Types"
#'   )
#'
#'   # --- Example 3: Interactive type selection ---
#'   # This block will only run if the session is interactive.
#'   if (interactive()) {
#'     message("Running interactive example. You will be prompted for selection.")
#'     fig_interactive <- plotByType(data = sampleData)
#'   }
#'
#'   # In a real session, you would display the plot by typing its name, e.g., `fig1`.
#' }
plotByType <- function(data, interactionTypes = NULL, showLabel = FALSE,
                       nodeColor = NULL, title = NULL) {

  # --- 1. Validate inputs ---
  validateTypePlotInputs(data)

  # --- 2. Handle interactive type selection ---
  if (is.null(interactionTypes)) {
    # This helper returns all types if non-interactive
    selectedTypes <- selectTypesInteractively(data)
    if (is.null(selectedTypes)) {
      warning("No interaction types selected. Returning NULL.")
      return(invisible(NULL))
    }
  } else {
    selectedTypes <- interactionTypes
  }

  # --- 3. Prepare data based on selection ---
  filteredData <- data[data$type %in% selectedTypes, ]
  if (nrow(filteredData) == 0) {
    warning("No data remains after filtering for specified interaction types.")
    return(invisible(NULL))
  }

  # --- 4. Prepare and create the Sankey plot ---
  sankeyData <- prepareTypePlotSankeyData(filteredData)
  sankeyData <- setTypePlotSankeyColors(sankeyData, nodeColor)
  createTypePlotSankey(sankeyData, showLabel, title)
}

# --- Helper Functions ---

#' Validate Inputs for Plot by Interaction Types
#' @noRd
validateTypePlotInputs <- function(data) {
  requiredCols <- c("intersection", "term_name", "type", "lncRNAId")
  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' must contain required columns" = all(requiredCols %in% colnames(data))
  )
}

#' Interactively Select Interaction Types or Default to All
#' @noRd
selectTypesInteractively <- function(data) {
  availableTypes <- unique(data$type)
  if (!interactive()) {
    return(availableTypes)
  }

  message("Select interaction types to plot:")
  choices <- c("All types", availableTypes)
  choicesText <- paste0(seq_along(choices), ". ", choices)
  message(paste(choicesText, collapse = "\n"))

  prompt <- "Enter numbers (comma-separated, e.g., 1 or 2,3): "
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
    return(availableTypes)
  }

  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(choices)]
  if (length(validIndices) == 0) {
    warning("Selected numbers do not correspond to available types.")
    return(NULL)
  }

  return(availableTypes[validIndices - 1])
}

#' Prepare Data Structure for Type Sankey Plot
#' @noRd
prepareTypePlotSankeyData <- function(data) {
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

#' Set Colors for Type Sankey Plot
#' @noRd
setTypePlotSankeyColors <- function(sankeyData, nodeColor) {
  node <- sankeyData$node
  link <- sankeyData$link
  if (is.null(nodeColor)) {
    nNodes <- nrow(node)
    palette <- Polychrome::createPalette(max(nNodes, 3), c("#ff0000", "#00ff00", "#0000ff"))
    node$color <- palette[seq_len(nNodes)]
  } else {
    node$color <- nodeColor
  }
  link$color <- node$color[link$IDtarget + 1]
  list(node = node, link = link)
}

#' Create Type Sankey Plot with Plotly
#' @noRd
createTypePlotSankey <- function(sankeyData, showLabel, title) {
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
