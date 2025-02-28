#' Plot Sankey Diagram by Interaction Type
#'
#' Creates a Sankey diagram visualizing interactions from a gProfiler output table, 
#' filtered by interaction type (cis/trans). Allows customization of labels, colors, and titles.
#'
#' @param data data.frame. A data frame containing gProfiler interaction results.
#'        Must include columns: `term_name`, `intersection`, `lncRNA.id`, `type`, and optionally `source`.
#' @param cis logical, optional. If `TRUE` (default), includes cis interactions.
#' @param trans logical, optional. If `TRUE` (default), includes trans interactions.
#' @param label logical, optional. If `TRUE`, shows node labels. If `FALSE` (default), no labels.
#' @param color character, optional. Single color for all nodes and links. If `NULL` (default), uses a random palette.
#' @param title character, optional. Diagram title. If `NULL` (default), no title.
#'
#' @return plotly object. A Sankey diagram visualizing filtered interactions.
#'
#' @examples
#' # Assuming 'combined_table' is your gProfiler interaction data frame
#' # fig <- plot_by_action(data = combined_table, cis = TRUE, trans = FALSE)
#' # fig # Display the plot
#'
#' @import plotly
#' @import Polychrome
#' @import dplyr
#' @export
plot_by_action <- function(data, cis = TRUE, trans = TRUE, label = FALSE, color = NULL, title = NULL) {
  # Validate inputs and filter data
  data <- validate_inputs_action_plot(data, cis, trans)
  
  # Prepare Sankey data
  sankey_data <- prepare_data_action_plot(data)
  
  # Set colors for nodes and links
  sankey_data <- set_sankey_colors_action_plot(sankey_data, color)
  
  # Create the Sankey plot
  fig <- create_sankey_action_plot(sankey_data, label, title)
  
  return(fig)
}

# --- Helper Functions  ---
#' Validate inputs for plot by interactions
#'
#' Ensures the input data has required columns and filters interactions based on type.
#'
#' @param data data.frame. Input data frame.
#' @param cis logical. Include cis interactions.
#' @param trans logical. Include trans interactions.
#' @return data.frame. Filtered data.
#' @noRd
validate_inputs_action_plot <- function(data, cis, trans) {
  required_columns <- c("intersection", "term_name", "type", "lncRNA.id")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(required_columns, collapse = ", "))
  }
  if (!cis && !trans) {
    stop("At least one of 'cis' or 'trans' must be TRUE.")
  }
  if (!cis) {
    data <- data[data$type != "cis", ]
  }
  if (!trans) {
    data <- data[data$type != "trans", ]
  }
  if (nrow(data) == 0) {
    stop("No data to plot after filtering.")
  }
  return(data)
}

#' Prepare Sankey Data
#'
#' Constructs link and node data frames from filtered interaction data.
#'
#' @param data data.frame. Filtered data frame.
#' @return list. Contains `link` and `node` data frames for Sankey plotting.
#' @noRd
prepare_data_action_plot <- function(data) {
  GOprot <- data.frame(source = data$intersection, target = data$term_name, stringsAsFactors = FALSE)
  LNCprot <- data.frame(source = data$lncRNA.id, target = data$intersection, stringsAsFactors = FALSE)
  
  link <- rbind(GOprot, LNCprot)
  link$value <- 1
  
  node <- data.frame(name = unique(c(link$source, link$target)), stringsAsFactors = FALSE)
  link$IDsource <- match(link$source, node$name) - 1
  link$IDtarget <- match(link$target, node$name) - 1
  list(link = link, node = node)
}

#' Set Sankey Colors
#'
#' Assigns colors to nodes and links in the Sankey diagram.
#'
#' @param sankey_data list. Contains `link` and `node` data frames.
#' @param color character. User-specified color or NULL for random palette.
#' @return list. Updated `sankey_data` with color assignments.
#' @noRd
set_sankey_colors_action_plot <- function(sankey_data, color) {
  node <- sankey_data$node
  link <- sankey_data$link
  if (is.null(color)) {
    node$color <- Polychrome::createPalette(nrow(node), c("#ff0000", "#00ff00", "#0000ff"))
  } else {
    node$color <- color
  }
  link$color <- node$color[link$IDtarget + 1]
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
create_sankey_action_plot <- function(sankey_data, label, title) {
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
