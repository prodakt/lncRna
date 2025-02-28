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
  validate_inputs(data, select_terms, color_selected, by = "terms")
  if (is.null(select_terms)) {
    select_terms <- interactive_selection(data, by = "terms")
  }
  prepared_data <- prepare_data(data, select_terms, color_selected, by = "terms")
  sankey_data <- prepare_sankey_data(prepared_data$data)
  sankey_data <- set_sankey_colors(sankey_data, color, color_selected, prepared_data$pathway, by = "terms")
  fig <- create_sankey_plot(sankey_data, label, title)
  return(fig)
}

#' Plot Sankey Diagram for gProfiler Target Connections
#'
#' Creates a Sankey diagram visualizing the connections of nodes to selected target protein-coding genes from a gProfiler output table.
#'
#' @param data data.frame. A data frame containing combined output from trans/cis_interactions functions.
#'        Must include columns: `term_name`, `intersection`, `lncRNA.id`, `type`, and optionally `source`.
#' @param select_target character vector, optional. Target gene IDs from `intersection` to highlight or plot. If `NULL` (default), prompts interactive selection.
#' @param label logical, optional. If `TRUE`, shows node labels. If `FALSE` (default), no labels.
#' @param color character, optional. Single color for all nodes/links. If `NULL` (default), uses a random palette.
#' @param title character, optional. Diagram title. If `NULL` (default), no title.
#' @param color_selected logical, optional. With `select_target`, if `TRUE`, plots all data but colors only selected targets. If `FALSE` (default), plots only selected targets.
#'
#' @return plotly object. Sankey diagram of target connections.
#' @export
plot_by_target <- function(data, select_target = NULL, label = FALSE, color = NULL, title = NULL, color_selected = FALSE) {
  validate_inputs(data, select_target, color_selected, by = "target")
  if (is.null(select_target)) {
    select_target <- interactive_selection(data, by = "target")
  }
  prepared_data <- prepare_data(data, select_target, color_selected, by = "target")
  sankey_data <- prepare_sankey_data(prepared_data$data)
  sankey_data <- set_sankey_colors(sankey_data, color, color_selected, prepared_data$pathway, by = "target")
  fig <- create_sankey_plot(sankey_data, label, title)
  return(fig)
}

# --- Helper Functions  ---

#' Validate Inputs for Plot Functions
#'
#' @param data data.frame. Input data frame.
#' @param select_items character vector. Selected items (terms or targets).
#' @param color_selected logical. Color selection flag.
#' @param by character. Selection type: "terms" or "targets".
#' @noRd
validate_inputs <- function(data, select_items, color_selected, by = "terms") {
  required_columns <- c("intersection", "term_name", "type", "lncRNA.id")
  if (!all(required_columns %in% colnames(data))) {
    stop("Input 'data' must contain columns: ", paste(required_columns, collapse = ", "))
  }
  if (color_selected && is.null(select_items)) {
    stop("To use 'color_selected = TRUE', please provide 'select_", by, "'.")
  }
  if (!is.null(select_items)) {
    column_to_check <- if (by == "terms") "term_name" else "intersection"
    if (!all(select_items %in% data[[column_to_check]])) {
      stop("Selected ", by, " not present in provided data.")
    }
  }
  if (nrow(data) == 0) {
    stop("No data to plot.")
  }
}

#' Interactive Selection for Terms or Targets
#'
#' @param data data.frame. Input data frame.
#' @param by character. Selection type: "terms" or "targets".
#' @return character vector. Selected items.
#' @noRd
interactive_selection <- function(data, by = "terms") {
  column_to_select <- if (by == "terms") "term_name" else "intersection"
  available_items <- unique(data[[column_to_select]])
  
  cat("Select ", by, " to plot:\n")
  cat("1. All ", by, "\n")
  for (i in seq_along(available_items)) {
    cat(paste0(i + 1, ". ", available_items[i], "\n"))
  }
  
  selected_indices_str <- readline(prompt = "Enter the number(s) (comma-separated, e.g., 1 or 2,3): ")
  selected_indices <- as.integer(strsplit(selected_indices_str, ",")[[1]])
  
  if (length(selected_indices) == 0 || any(is.na(selected_indices))) {
    stop("Invalid input: please enter comma-separated numbers.")
  }
  
  if (1 %in% selected_indices) {
    cat("You selected: All ", by, "\n")
    return(available_items)
  } else {
    item_indices <- selected_indices - 1
    if (any(item_indices < 1) || any(item_indices > length(available_items))) {
      stop("Invalid selection: please enter numbers between 1 and ", length(available_items) + 1, ".")
    }
    select_items <- available_items[item_indices]
    cat("You selected ", by, ": ", paste(select_items, collapse = ", "), "\n")
    return(select_items)
  }
}

#' Prepare Data for Sankey Diagram
#'
#' @param data data.frame. Input data frame.
#' @param select_items character vector. Selected items (terms or targets).
#' @param color_selected logical. Color selection flag.
#' @param by character. Selection type: "terms" or "targets".
#' @return list. Contains 'data' (filtered or full) and 'pathway' (selected data if color_selected is TRUE).
#' @noRd
prepare_data <- function(data, select_items, color_selected, by = "terms") {
  column_to_select <- if (by == "terms") "term_name" else "intersection"
  if (color_selected) {
    pathway <- data[data[[column_to_select]] %in% select_items, ]
    return(list(data = data, pathway = pathway))
  } else {
    filtered_data <- data[data[[column_to_select]] %in% select_items, ]
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
#' @param sankey_data list. Sankey data with link and node data frames.
#' @param color character. User-specified color.
#' @param color_selected logical. Color selection flag.
#' @param pathway data.frame. Selected data if color_selected is TRUE.
#' @param by character. Selection type: "terms" or "targets".
#' @return list. Updated Sankey data with colors.
#' @noRd
set_sankey_colors <- function(sankey_data, color, color_selected, pathway, by = "terms") {
  node <- sankey_data$node
  link <- sankey_data$link
  
  if (is.null(color)) {
    node$color <- Polychrome::createPalette(nrow(node), c("#ff0000", "#00ff00", "#0000ff"))
  } else {
    node$color <- color
  }
  
  if (color_selected && !is.null(pathway)) {
    pathway_nodes <- unique(c(pathway$intersection, pathway$term_name, pathway$lncRNA.id))
    node$color[!node$name %in% pathway_nodes] <- "gray"
  }
  
  link$color <- node$color[link$IDtarget + 1]
  
  if (color_selected && !is.null(pathway)) {
    column_to_check <- if (by == "terms") "term_name" else "intersection"
    link$color[!link$target %in% pathway[[column_to_check]]] <- "gray"
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
  
  return(fig)
}
