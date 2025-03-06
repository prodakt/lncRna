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
plot_by_type <- function(data, type, label = FALSE, color = NULL, title = NULL) {
  # Validate inputs and filter data
  type <- validate_inputs_type_plot(data, type)
  
  # Prepare Sankey data
  sankey_data <- prepare_data_type_plot(data)
  
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
  
  # Handle missing type parameter - interactive selection
  if (missing(type) || is.null(type)) {
    unique_types <- unique(data$type)
    if(length(unique_types) == 0) {
      stop("No interaction types found in 'type' column of the data.")
    }
    cat("Available interaction types in 'type' column:\n")
    for (i in seq_along(unique_types)) {
      cat(paste0(i, ". ", unique_types[i], "\n"))
    }
    type_choice_index <- readline(prompt = "Choose type number from console or type 'all' for all types: ")
    
    if(tolower(type_choice_index) == 'all'){
      type <- unique_types
    } else {
      type_choice_index <- as.integer(type_choice_index)
      if (is.na(type_choice_index) || type_choice_index < 1 || type_choice_index > length(unique_types)) {
        stop("Invalid input. Please choose a number from the list or type 'all'.")
      }
      type <- unique_types[type_choice_index]
    }
    cat("You selected interaction type(s): ", paste(type, collapse = ", "), "\n")
  } else {
    # Validate type parameter if provided
    if (!is.character(type)) {
      stop("'type' must be a character vector specifying interaction types.")
    }
    if (length(type) == 0) {
      stop("'type' must be a non-empty character vector specifying interaction types.")
    }
    # Check if all specified type types exist in data
    if (!all(type %in% unique(data$type))) {
      missing_types <- type[!type %in% unique(data$type)]
      stop("Specified interaction types not found in data: ", paste(missing_types, collapse = ", "))
    }
  }
  
  
  # Filter data based on type types
  data <- data[data$type %in% type, ]
  
  if (nrow(data) == 0) {
    stop("No data to plot after filtering with specified types.")
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
prepare_data_type_plot <- function(data) {
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
