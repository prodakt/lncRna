#' Plot Sankey Diagram for gProfiler Term Connections
#'
#' Creates a Sankey diagram visualizing the connections of nodes to selected terms from a gProfiler output table.
#' This function is designed to display relationships between lncRNAs, intersection (e.g., genomic regions),
#' and functional terms (like GO terms or pathways) obtained from gProfiler analyses.
#'
#' @param data data.frame. A data frame containing combined output from trans/cis_interactions functions.
#'        This data frame must include at least the following columns:
#'        `term_name` (functional term name from gProfiler),
#'        `intersection` (e.g., genomic region or gene intersection),
#'        `lncRNA.id` (lncRNA identifier),
#'        `type` (interaction type, e.g., "trans", "cis"),
#'        and optionally `source` (data source identifier).
#' @param select_terms character vector, optional. A vector of term names from the `term_name` column of the input `data`
#'        to be highlighted or exclusively plotted. If `NULL` (default), the function will prompt the user to interactively select terms.
#' @param label logical, optional. If `TRUE`, node labels (term names, intersection IDs, lncRNA IDs) are displayed on the Sankey diagram nodes.
#'        If `FALSE` (default), no labels are shown.
#' @param color character, optional. A single color name (e.g., "blue", "steelblue") to manually apply to all nodes and links in the diagram.
#'        If `NULL` (default), a random color palette is generated for nodes, and links inherit colors from their target nodes.
#' @param title character, optional. The title of the Sankey diagram. If `NULL` (default), no title is displayed.
#' @param color_selected logical, optional. Used in conjunction with `select_terms`.
#'        If `TRUE`, all interactions from the input `data` are plotted, but only nodes and links associated with the `select_terms`
#'        are colored according to the `color` parameter (or default palette). Nodes and links not associated with `select_terms` are colored gray.
#'        If `FALSE` (default), only interactions directly involving the `select_terms` are plotted, and coloring applies to these selected elements.
#'        Cannot be `TRUE` if `select_terms = NULL`.
#'
#' @return plotly object.
#' Returns a plotly Sankey diagram object visualizing term connections.
#' This object can be further customized using plotly functions or directly rendered.
#'
#' @examples
#' # Assuming 'combined_table' is your data frame containing gProfiler interaction results
#' # (replace 'combined_table' with your actual data frame name)
#'
#' # Example 1: Basic Sankey diagram plotting all terms with interactive selection
#' # fig1 <- plot_by_terms(data = combined_table)
#' # fig1 # Display the plot
#'
#' # Example 2: Plotting only interactions for specific terms (e.g., "response to stress")
#' # fig2 <- plot_by_terms(data = combined_table, select_terms = "response to stress")
#' # fig2 # Display the plot
#'
#' # Example 3: Plotting interactions for multiple selected terms
#' # fig3 <- plot_by_terms(data = combined_table, select_terms = c("response to stress", "apoptotic process"))
#' # fig3 # Display the plot
#'
#' # Example 4: Sankey diagram with node labels enabled
#' # fig4 <- plot_by_terms(data = combined_table, label = TRUE)
#' # fig4 # Display the plot
#'
#' # Example 5: Sankey diagram with a custom title
#' # fig5 <- plot_by_terms(data = combined_table, title = "Sankey Diagram of Term Interactions")
#' # fig5 # Display the plot
#'
#' # Example 6: Sankey diagram with a single color for all elements (e.g., "skyblue")
#' # fig6 <- plot_by_terms(data = combined_table, color = "skyblue")
#' # fig6 # Display the plot
#'
#' # Example 7: Highlighting interactions related to "response to stress" term, while showing all data
#' # fig7 <- plot_by_terms(data = combined_table, select_terms = "response to stress", color_selected = TRUE)
#' # fig7 # Display the plot
#'
#' # Example 8: Highlighting with a custom color for selected terms
#' # fig8 <- plot_by_terms(data = combined_table, select_terms = "response to stress", color_selected = TRUE, color = "forestgreen")
#' # fig8 # Display the plot
#'
#' # Example 9: Interactive term selection from console
#' # fig9 <- plot_by_terms(data = combined_table, select_terms = NULL)
#' # fig9 # Display the plot
#'
#' @import plotly
#' @import Polychrome
#' @import dplyr
#' @export
plot_by_terms <- function(data, select_terms = NULL, label = FALSE, color = NULL, title = NULL, color_selected = FALSE) {
  # Validate inputs
  validate_inputs_plot_terms(data, select_terms, color_selected)
  
  # Handle interactive term selection if select_terms is NULL
  if (is.null(select_terms)) {
    select_terms <- interactive_term_selection(data)
  }
  
  # Prepare data based on selection
  prepared_data <- prepare_data(data, select_terms, color_selected)
  
  # Prepare Sankey data
  sankey_data <- prepare_sankey_data(prepared_data$data)
  
  # Set colors for nodes and links
  sankey_data <- set_sankey_colors(sankey_data, color, color_selected, prepared_data$pathway)
  
  # Create the Sankey plot
  fig <- create_sankey_plot(sankey_data, label, title)
  
  return(fig)
}

# --- Helper Functions  ---

#' Validate Inputs for plot_by_terms
#'
#' @param data data.frame. Input data frame.
#' @param select_terms character vector. Selected terms.
#' @param color_selected logical. Color selection flag.
#' @noRd
validate_inputs_plot_terms <- function(data, select_terms, color_selected) {
  required_columns <- c( "intersection", "term_name", "type", "lncRNA.id")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(required_columns, collapse = ", "))
  }
  if (color_selected && is.null(select_terms)) {
    stop("To use 'color_selected = TRUE', please provide 'select_terms'.")
  }
  if (!is.null(select_terms) && !all(select_terms %in% data$term_name)) {
    stop("Selected terms not present in provided data.")
  }
  if (nrow(data) == 0) {
    stop("No data to plot.")
  }
}

#' Interactive Term Selection
#'
#' @param data data.frame. Input data frame.
#' @return character vector. Selected terms.
#' @noRd
interactive_term_selection <- function(data) {
  # Get unique terms from the term_name column
  available_terms <- unique(data$term_name)
  
  # Display selection options starting with "1. All terms"
  cat("Select terms to plot:\n")
  cat("1. All terms\n")
  for (i in seq_along(available_terms)) {
    cat(paste0(i + 1, ". ", available_terms[i], "\n"))
  }
  
  # Prompt user for input
  selected_indices_str <- readline(prompt = "Enter the number(s) (comma-separated, e.g., 1 or 2,3): ")
  selected_indices <- as.integer(strsplit(selected_indices_str, ",")[[1]])
  
  # Validate input
  if (length(selected_indices) == 0 || any(is.na(selected_indices))) {
    stop("Invalid input: please enter comma-separated numbers.")
  }
  
  # Check if "1" (All terms) is selected
  if (1 %in% selected_indices) {
    cat("You selected: All terms\n")
    return(available_terms)
  } else {
    # Adjust indices since individual terms start at 2
    term_indices <- selected_indices - 1
    if (any(term_indices < 1) || any(term_indices > length(available_terms))) {
      stop("Invalid selection: please enter numbers between 1 and ", length(available_terms) + 1, ".")
    }
    select_terms <- available_terms[term_indices]
    cat("You selected terms: ", paste(select_terms, collapse = ", "), "\n")
    return(select_terms)
  }
}

#' Prepare Data for Sankey Diagram
#'
#' @param data data.frame. Input data frame.
#' @param select_terms character vector. Selected terms.
#' @param color_selected logical. Color selection flag.
#' @return list. Contains 'data' (filtered or full data) and 'pathway' (selected data if color_selected is TRUE).
#' @noRd
prepare_data <- function(data, select_terms, color_selected) {
  if (color_selected) {
    pathway <- data[data$term_name %in% select_terms, ]
    return(list(data = data, pathway = pathway))
  } else {
    filtered_data <- data[data$term_name %in% select_terms, ]
    return(list(data = filtered_data, pathway = NULL))
  }
}

#' Prepare Sankey Data
#'
#' @param data data.frame. Prepared data frame.
#' @return list. Sankey data with link and node data frames.
#' @noRd
prepare_sankey_data <- function(data) {
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
#' @param sankey_data list. Sankey data with link and node data frames.
#' @param color character. User-specified color.
#' @param color_selected logical. Color selection flag.
#' @param pathway data.frame. Selected data if color_selected is TRUE.
#' @return list. Updated Sankey data with colors.
#' @noRd
set_sankey_colors <- function(sankey_data, color, color_selected, pathway) {
  node <- sankey_data$node
  link <- sankey_data$link
  
  # Assign node colors
  if (is.null(color)) {
    n_nodes <- nrow(node)
    # Generate palette and ensure exact length
    palette <- Polychrome::createPalette(max(n_nodes, 3), c("#ff0000", "#00ff00", "#0000ff"))
    print(length(palette))  # Debug output
    node$color <- palette[1:n_nodes]  # Take exactly n_nodes colors
  } else {
    node$color <- color
  }
  
  # Adjust colors if color_selected is TRUE
  if (color_selected && !is.null(pathway)) {
    pathway_nodes <- unique(c(pathway$intersection, pathway$term_name, pathway$lncRNA.id))
    node$color[!node$name %in% pathway_nodes] <- "gray"
  }
  
  # Assign link colors based on target nodes
  link$color <- node$color[link$IDtarget + 1]
  
  # If color_selected, gray out links not associated with selected terms
  if (color_selected && !is.null(pathway)) {
    link$color[!link$target %in% pathway$term_name] <- "gray"
  }
  
  list(node = node, link = link)
}

#' Create Sankey Plot
#'
#' @param sankey_data list. Sankey data with link and node data frames.
#' @param label logical. Whether to display labels.
#' @param title character. Plot title.
#' @return plotly object. Sankey plot.
#' @noRd
create_sankey_plot <- function(sankey_data, label, title) {
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
  
  fig
}
