#' Plot Sankey Diagram for gProfiler Term Connections
#'
#' Creates a Sankey diagram visualizing connections to selected gProfiler terms.
#' Allows for interactive term selection in an interactive session.
#'
#' @param data A data.frame from `trans/cis_interactions` functions. Must
#'   include columns: `term_name`, `intersection`, `lncRNAId`, and `type`.
#' @param selectedTerms A character vector of term names to highlight or plot.
#'   If `NULL` (default) and the session is interactive, the user will be
#'   prompted to select terms. In a non-interactive session, `NULL` results
#'   in all terms being plotted.
#' @param showLabel A logical value; if `TRUE`, node labels are displayed.
#'   Defaults to `FALSE`.
#' @param nodeColor A single color name (e.g., "blue") for all nodes. If `NULL`
#'   (default), a random palette is used.
#' @param title A character string for the diagram title. Defaults to `NULL`.
#' @param colorSelected A logical value. If `TRUE`, all interactions are
#'   plotted, but only those related to `selectedTerms` are colored (others are
#'   gray). If `FALSE` (default), only interactions involving `selectedTerms` are
#'   plotted. Requires `selectedTerms` to be specified.
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
#'   intersection = c("GeneA", "GeneB", "GeneC", "GeneA"),
#'   lncRNAId = c("Lnc1", "Lnc1", "Lnc2", "Lnc2"),
#'   type = c("cis", "trans", "cis", "trans"),
#'   stringsAsFactors = FALSE
#' )
#'
#' if (requireNamespace("plotly", quietly = TRUE) &&
#'     requireNamespace("Polychrome", quietly = TRUE)) {
#'
#'   # --- Example 1: Non-interactive plot for a specific term ---
#'   fig1 <- plotByTerms(
#'     data = sampleData,
#'     selectedTerms = "GO:2",
#'     title = "Interactions for GO:2"
#'   )
#'
#'   # --- Example 2: Non-interactive plot highlighting a specific term ---
#'   fig2 <- plotByTerms(
#'     data = sampleData,
#'     selectedTerms = "GO:1",
#'     colorSelected = TRUE,
#'     title = "All Interactions (GO:1 Highlighted)"
#'   )
#'
#'   # --- Example 3: Interactive term selection ---
#'   # This block will only run if the session is interactive.
#'   if (interactive()) {
#'     message("Running interactive example. You will be prompted for selection.")
#'     fig_interactive <- plotByTerms(data = sampleData)
#'   }
#'
#'   # In a real session, you would display the plot by typing its name, e.g., `fig1`.
#' }
plotByTerms <- function(data, selectedTerms = NULL, showLabel = FALSE,
                        nodeColor = NULL, title = NULL, colorSelected = FALSE) {

  # --- 1. Validate inputs ---
  validateTermPlotInputs(data, selectedTerms, colorSelected)

  # --- 2. Handle interactive term selection ---
  if (is.null(selectedTerms)) {
    # This helper returns all terms if non-interactive
    selectedTerms <- selectTermsInteractively(data)
    if (is.null(selectedTerms)) {
      warning("No terms selected. Returning NULL.")
      return(invisible(NULL))
    }
  }

  # --- 3. Prepare data based on selection ---
  preparedData <- prepareTermPlotData(data, selectedTerms, colorSelected)
  if (is.null(preparedData)) {
    return(invisible(NULL))
  }

  # --- 4. Prepare and create the Sankey plot ---
  sankeyData <- prepareTermPlotSankeyData(preparedData$data)
  sankeyData <- setTermPlotSankeyColors(
    sankeyData, nodeColor, colorSelected, preparedData$pathway
  )
  createTermPlotSankey(sankeyData, showLabel, title)
}

# --- Helper Functions ---

#' Validate Inputs for plotByTerms
#' @noRd
validateTermPlotInputs <- function(data, selectedTerms, colorSelected) {
  requiredCols <- c("intersection", "term_name", "type", "lncRNAId")
  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' must contain required columns" = all(requiredCols %in% colnames(data))
  )
  if (colorSelected && is.null(selectedTerms)) {
    stop("To use 'colorSelected = TRUE', please provide 'selectedTerms'.")
  }
  if (!is.null(selectedTerms) && !all(selectedTerms %in% data$term_name)) {
    # Check for any valid terms before stopping
    if (!any(selectedTerms %in% data$term_name)) {
      stop("None of the selected terms were found in the data.")
    } else {
      warning("Some selected terms were not found in the data and will be ignored.")
    }
  }
}

#' Interactively Select Terms or Default to All
#' @noRd
selectTermsInteractively <- function(data) {
  availableTerms <- unique(data$term_name)
  if (!interactive()) {
    return(availableTerms)
  }

  message("Select terms to plot:")
  choices <- c("All terms", availableTerms)
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
    return(availableTerms)
  }

  validIndices <- selectedIndices[selectedIndices > 1 & selectedIndices <= length(choices)]
  if (length(validIndices) == 0) {
    warning("Selected numbers do not correspond to available terms.")
    return(NULL)
  }

  return(availableTerms[validIndices - 1])
}

#' Prepare Data for Sankey Diagram
#' @noRd
prepareTermPlotData <- function(data, selectedTerms, colorSelected) {
  # Filter out any selected terms not in data
  validSelectedTerms <- intersect(selectedTerms, data$term_name)
  if (length(validSelectedTerms) == 0) {
    warning("None of the selected terms are in the data after filtering.")
    return(NULL)
  }

  if (colorSelected) {
    pathway <- data[data$term_name %in% validSelectedTerms, ]
    return(list(data = data, pathway = pathway))
  } else {
    filteredData <- data[data$term_name %in% validSelectedTerms, ]
    return(list(data = filteredData, pathway = NULL))
  }
}

#' Prepare Data Structure for Sankey Plot
#' @noRd
prepareTermPlotSankeyData <- function(data) {
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
setTermPlotSankeyColors <- function(sankeyData, nodeColor, colorSelected, pathway) {
  node <- sankeyData$node
  link <- sankeyData$link

  if (is.null(nodeColor)) {
    nNodes <- nrow(node)
    palette <- Polychrome::createPalette(max(nNodes, 3), c("#ff0000", "#00ff00", "#0000ff"))
    node$color <- palette[seq_len(nNodes)]
  } else {
    node$color <- nodeColor
  }

  if (colorSelected && !is.null(pathway)) {
    pathwayNodes <- unique(c(pathway$intersection, pathway$term_name, pathway$lncRNAId))
    node$color[!node$name %in% pathwayNodes] <- "gray"
  }

  link$color <- node$color[link$IDtarget + 1]

  if (colorSelected && !is.null(pathway)) {
    link$color[!link$target %in% pathway$term_name] <- "gray"
  }

  list(node = node, link = link)
}

#' Create Sankey Plot with Plotly
#' @noRd
createTermPlotSankey <- function(sankeyData, showLabel, title) {
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
