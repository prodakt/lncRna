#' Plot Sankey Diagram by Specified Interaction Types
#'
#' Creates a Sankey diagram visualizing interactions from a gProfiler output table,
#' filtered by user-specified interaction types, with interactive type selection from the console if not specified.
#' Allows customization of labels, colors, and titles.
#'
#' @param data data.frame. A data frame containing gProfiler interaction results.
#'        Must include columns: `term_name`, `intersection`, `lncRNA.id`, `type`, and optionally `source`.
#' @param type character vector, optional. Specifies the interaction types to include in the plot (e.g., "cis", "trans", "TAR").
#'        Must match values in the `type` column of `data`. If missing, user will be prompted to choose from console.
#' @param label logical, optional. If `TRUE`, shows node labels. If `FALSE` (default), no labels are shown.
#' @param color character, optional. Single color for all nodes and links. If `NULL` (default), uses a random palette.
#' @param title character, optional. Diagram title. If `NULL` (default), no title is displayed.
#'
#' @return plotly object. A Sankey diagram visualizing filtered interactions.
#'
#' @examples
#' # Assuming 'combined_table' is your gProfiler interaction data frame
#'
#' # Example 1: Plot only "cis" interactions
#' # fig1 <- plot_by_type(data = combined_table, type = "cis")
#' # fig1 # Display the plot
#'
#' # Example 2: Plot "cis" and "trans" interactions
#' # fig2 <- plot_by_type(data = combined_table, type = c("cis", "trans"), label = TRUE)
#' # fig2 # Display the plot
#'
#' # Example 3: Plot custom interaction type with a title and color
#' # fig3 <- plot_by_type(data = combined_table, type = "TAR", title = "TAR Interactions", color = "blue")
#' # fig3 # Display the plot
#'
#' # Example 4: Interactive type selection from console
#' # fig4 <- plot_by_type(data = combined_table) # User will be prompted to select type
#' # fig4 # Display the plot
#'
#' @import plotly
#' @import Polychrome
#' @import dplyr
#' @export
plot_by_type <- function(data, type = NULL, label = FALSE, color = NULL, title = NULL, color_selected = FALSE) {
  # Validate inputs and filter data
  validate_inputs_type_plot(data, type)
  
  # Handle interactive term selection if select_terms is NULL
  if (is.null(type)) {
    selected_type <- interactive_types_selection(data)
  }
  
  # Prepare data based on selection
  prepared_data <- prepare_data_type(data, selected_type, color_selected)
  
  # Prepare Sankey data
  sankey_data <- prepare_data_type_plot(prepared_data$data)
  
  # Set colors for nodes and links
  sankey_data <- set_sankey_colors_type_plot(sankey_data, color)
  
  # Create the Sankey plot
  fig <- create_sankey_type_plot(sankey_data, label, title)
  
  return(fig)
}

# --- Helper Functions ---

#' Validate Inputs for Plot by Interaction Types
#'
#' Ensures the input data has required columns and filters interactions based on user-specified types.
#' If type is not specified, prompts user to select from console.
#'
#' @param data data.frame. Input data frame.
#' @param type character vector. Interaction types to include.
#' @return data.frame. Filtered data.
#' @noRd
validate_inputs_type_plot <- function(data, type) {
  # Validate required columns
  required_columns <- c("intersection", "term_name", "type", "lncRNA.id")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(required_columns, collapse = ", "))
  }
  
}

#' Interactive Type Selection
#'
#' @param data data.frame. Input data frame.
#' @return character vector. Selected types.
#' @noRd
interactive_types_selection <- function(data) {
  # Get unique types from the type column
  available_types <- unique(data$type)
  
  # Display selection options starting with "1. All types"
  cat("Select type to plot:\n")
  cat("1. All types\n")
  for (i in seq_along(available_types)) {
    cat(paste0(i + 1, ". ", available_types[i], "\n"))
  }
  
  # Prompt user for input
  selected_indices_str <- readline(prompt = "Enter the number(s) (comma-separated, e.g., 1 or 2,3): ")
  selected_indices <- as.integer(strsplit(selected_indices_str, ",")[[1]])
  
  # Validate input
  if (length(selected_indices) == 0 || any(is.na(selected_indices))) {
    stop("Invalid input: please enter comma-separated numbers.")
  }
  
  # Check if "1" (All types) is selected
  if (1 %in% selected_indices) {
    cat("You selected: All types\n")
    return(available_types)
  } else {
    # Adjust indices since individual types start at 2
    type_indices <- selected_indices - 1
    if (any(type_indices < 1) || any(type_indices > length(available_types))) {
      stop("Invalid selection: please enter numbers between 1 and ", length(available_types) + 1, ".")
    }
    select_types <- available_types[type_indices]
    cat("You selected types: ", paste(select_types, collapse = ", "), "\n")
    return(select_types)
  }
}

#' Prepare filtered Data for Sankey Diagram
#'
#' @param data data.frame. Input data frame.
#' @param select_type character vector. Selected types.
#' @param color_selected logical. Color selection flag.
#' @return list. Contains 'data' (filtered or full data) and 'pathway' (selected data if color_selected is TRUE).
#' @noRd
prepare_data_type <- function(data, select_type, color_selected) {
  if (color_selected) {
    pathway <- data[data$type %in% select_type, ]
    return(list(data = data, pathway = pathway))
  } else {
    filtered_data <- data[data$type %in% select_type, ]
    return(list(data = filtered_data, pathway = NULL))
  }
}

#' Prepare Sankey Data
#'
#' Constructs link and node data frames from filtered interaction data.
#'
#' @param sankey_data data.frame. Filtered data frame.
#' @return list. Contains `link` and `node` data frames for Sankey plotting.
#' @noRd
prepare_data_type_plot <- function(sankey_data) {
  GOprot <- data.frame(source = sankey_data$intersection, target = sankey_data$term_name, stringsAsFactors = FALSE)
  LNCprot <- data.frame(source = sankey_data$lncRNA.id, target = sankey_data$intersection, stringsAsFactors = FALSE)
  
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
set_sankey_colors_type_plot <- function(sankey_data, color) {
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
create_sankey_type_plot <- function(sankey_data, label, title) {
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
