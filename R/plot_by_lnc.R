#' Plot Sankey Diagram for Selected lncRNAs
#'
#' Creates a Sankey diagram visualizing interactions connected to specific lncRNAs from a gProfiler output table.
#' Allows highlighting or exclusive plotting of selected lncRNAs with customizable labeling, coloring, and titles.
#'
#' @param data data.frame. A data frame containing gProfiler interaction results.
#'        Must include columns: `term_name`, `intersection`, `lncRNA.id`, `type`, and optionally `source`.
#' @param select_lnc character vector, optional. Vector of lncRNA IDs to highlight or exclusively plot.
#'        If `NULL` (default), plots all interactions.
#' @param label logical, optional. If `TRUE`, shows node labels. If `FALSE` (default), no labels.
#' @param color character, optional. Single color for all nodes and links. If `NULL` (default), uses a random palette.
#' @param title character, optional. Diagram title. If `NULL` (default), no title.
#' @param color_selected logical, optional. If `TRUE`, plots all interactions but colors only those tied to `select_lnc`.
#'        If `FALSE` (default), plots only selected interactions. Requires `select_lnc`.
#'
#' @return plotly object. A Sankey diagram visualizing lncRNA interactions.
#'
#' @examples
#' # Assuming 'combined_table' is your gProfiler interaction data frame
#' # fig <- plot_by_lnc(data = combined_table, select_lnc = c("lncRNA1", "lncRNA2"), color_selected = TRUE)
#' # fig # Display the plot
#'
#' @import plotly
#' @import Polychrome
#' @import dplyr
#' @export
plot_by_lnc <- function(data, select_lnc = NULL, label = FALSE, color = NULL, title = NULL, color_selected = FALSE) {
  # Validate inputs and prepare data
  prepared_data <- validate_inputs_lnc_plot(data, select_lnc, color_selected)

  # Prepare Sankey data
  sankey_data <- prepare_data_lnc_plot(prepared_data$data)

  # Set colors for nodes and links
  sankey_data <- set_sankey_colors_lnc_plot(sankey_data, color, color_selected, prepared_data$pathway)

  # Create the Sankey plot
  fig <- create_sankey_lnc_plot(sankey_data, label, title)

  return(fig)
}

#' Validate and Prepare Data for lncRNA Plotting
#'
#' Ensures input data is valid and prepares it based on lncRNA selection and coloring options.
#'
#' @param data data.frame. Input data frame.
#' @param select_lnc character vector. Selected lncRNA IDs.
#' @param color_selected logical. Flag to determine coloring behavior.
#' @return list. Contains `data` (filtered or full) and `pathway` (selected data if applicable).
#' @noRd
validate_inputs_lnc_plot <- function(data, select_lnc, color_selected) {
  required_columns <- c("intersection", "term_name", "type", "lncRNA.id")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(required_columns, collapse = ", "))
  }
  if (color_selected && is.null(select_lnc)) {
    stop("To use 'color_selected = TRUE', please provide 'select_lnc'.")
  }
  if (!is.null(select_lnc)) {
    if (color_selected) {
      pathway <- data[data$lncRNA.id %in% select_lnc, ]
      if (nrow(pathway) == 0) {
        stop("Selected lncRNAs are not present in provided data.")
      }
      return(list(data = data, pathway = pathway))
    } else {
      filtered_data <- data[data$lncRNA.id %in% select_lnc, ]
      if (nrow(filtered_data) == 0) {
        stop("Selected lncRNAs are not present in provided data.")
      }
      return(list(data = filtered_data, pathway = NULL))
    }
  } else {
    return(list(data = data, pathway = NULL))
  }
}

#' Prepare Sankey Data
#'
#' Constructs link and node data frames from prepared interaction data.
#'
#' @param data data.frame. Prepared data frame.
#' @return list. Contains `link` and `node` data frames for Sankey plotting.
#' @noRd
prepare_data_lnc_plot <- function(data) {
  GOprot <- data.frame(source = data$intersection, target = data$term_name, stringsAsFactors = FALSE)
  LNCprot <- data.frame(source = data$lncRNA.id, target = data$intersection, stringsAsFactors = FALSE)
  link <- rbind(GOprot, LNCprot)
  link$value <- 1
  node <- data.frame(name = unique(c(link$source, link$target)), stringsAsFactors = FALSE)
  link$IDsource <- match(link$source, node$name) - 1
  link$IDtarget <- match(link$target, node$name) - 1
  list(link = link, node = node)
}

#' Set Sankey Colors for lncRNA Plotting
#'
#' Assigns colors to nodes and links, with optional graying out of unselected elements.
#'
#' @param sankey_data list. Contains `link` and `node` data frames.
#' @param color character. User-specified color or NULL for random palette.
#' @param color_selected logical. Flag to gray out unselected elements.
#' @param pathway data.frame. Data for selected lncRNAs if applicable.
#' @return list. Updated `sankey_data` with color assignments.
#' @noRd
set_sankey_colors_lnc_plot <- function(sankey_data, color, color_selected, pathway) {
  node <- sankey_data$node
  link <- sankey_data$link
  if (is.null(color)) {
    node$color <- Polychrome::createPalette(nrow(node), c("#ff0000", "#00ff00", "#0000ff"))
  } else {
    node$color <- color
  }
  if (color_selected && !is.null(pathway)) {
    selected_nodes <- unique(c(pathway$intersection, pathway$term_name, pathway$lncRNA.id))
    node$color[!node$name %in% selected_nodes] <- "gray"
  }
  link$color <- node$color[link$IDtarget + 1]
  if (color_selected && !is.null(pathway)) {
    selected_lnc <- pathway$lncRNA.id
    link$color[!link$source %in% selected_lnc & !link$target %in% selected_lnc] <- "gray"
  }
  list(node = node, link = link)
}

#' Create Sankey Plot
#'
#' Builds a Sankey diagram from prepared data using Plotly.
#'
#' @param sankey_data list. Contains `link` and `node` data frames with colors.
#' @param label logical. Whether to display node labels.
#' @param title character. Plot title.
#' @return plotly object. The Sankey diagram.
#' @noRd
create_sankey_lnc_plot <- function(sankey_data, label, title) {
  node <- sankey_data$node
  link <- sankey_data$link
  plot_node <- list(
    label = if (label) node$name else NULL,
    color = node$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5)
  )
  plot_link <- list(
    source = link$IDsource,
    target = link$IDtarget,
    value = link$value,
    color = link$color
  )
  fig <- plotly::plot_ly(
    type = "sankey",
    orientation = "h",
    node = plot_node,
    link = plot_link
  ) %>%
    plotly::layout(
      title = title,
      font = list(size = 10)
    )
  return(fig)
}
