#' Plot Functional Interactions as a Sankey Diagram
#'
#' Visualizes genomic interaction data (lncRNAId -> Target -> Functional Term)
#' as an interactive Sankey diagram using Plotly. Replaces `plot_by_terms`,
#' `plot_by_lnc`, etc., providing a unified interface for filtering and plotting.
#'
#' @param interactionData A `data.frame` containing interaction data, typically
#'   created by `annotateInteractions()`. Must contain columns: `"lncRNAId"`,
#'   `"intersection"`, `"term_name"`, and `"type"`.
#' @param groupBy A character string specifying the filtering criteria. Valid
#'   options are: `"lncRNAId"`, `"target"`, `"term"`, `"type"`. If `NULL`,
#'   interactive selection is triggered.
#' @param selectIds Optional character vector of IDs/terms to filter or highlight.
#'   If `NULL` and interactive, prompts for selection.
#' @param showLabels Logical. If `TRUE`, displays labels on nodes. Defaults to `FALSE`.
#' @param highlightSelected Logical. If `TRUE`, plots all data but colors only
#'   the nodes/links associated with `selectIds`. If `FALSE` (default), plots
#'   only the data matching `selectIds`.
#' @param color Optional character string. A single color for all nodes/links.
#'   If `NULL`, a random palette (Polychrome) is generated.
#' @param title Optional character string for the plot title.
#'
#' @return A `plotly` object containing the Sankey diagram.
#' @export
#' @importFrom plotly plot_ly layout
#' @importFrom Polychrome createPalette
#' @importFrom grDevices rgb
#'
#' @examples
#' # --- 1. Create mock interaction data ---
#' mockData <- data.frame(
#'   lncRNAId = c(rep("LNC1", 3), rep("LNC2", 2)),
#'   intersection = c("T1", "T2", "T1", "T3", "T2"),
#'   term_name = c("Stress", "Growth", "Stress", "Immunity", "Growth"),
#'   type = c("cis", "cis", "trans", "trans", "cis"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Run in non-interactive mode ---
#'
#' # Example 1: Filter by specific Term
#' fig1 <- plotSankeyInteractions(
#'   interactionData = mockData,
#'   groupBy = "term",
#'   selectIds = "Stress",
#'   title = "Stress Interactions"
#' )
#' # fig1
#'
#' # Example 2: Highlight specific lncRNAId
#' fig2 <- plotSankeyInteractions(
#'   interactionData = mockData,
#'   groupBy = "lncRNAId",
#'   selectIds = "LNC1",
#'   highlightSelected = TRUE,
#'   title = "Highlighting LNC1"
#' )
#' # fig2
#'
plotSankeyInteractions <- function(interactionData,
                                   groupBy = NULL,
                                   selectIds = NULL,
                                   showLabels = FALSE,
                                   highlightSelected = FALSE,
                                   color = NULL,
                                   title = NULL) {

  # --- 1. Input Validation ---
  interactionData <- as.data.frame(interactionData)

  required_cols <- c("lncRNAId", "intersection", "term_name", "type")
  if (!all(required_cols %in% colnames(interactionData))) {
    stop("Input 'interactionData' must contain columns: ",
         paste(required_cols, collapse = ", "))
  }

  group_map <- c(
    "lncRNAId" = "lncRNAId",
    "target" = "intersection",
    "term" = "term_name",
    "type" = "type"
  )

  # --- 2. Handle 'groupBy' Selection ---
  if (is.null(groupBy)) {
    if (interactive()) {
      groupBy <- .selectGroupByInteractive()
      if (is.null(groupBy)) return(invisible(NULL))
    } else {
      warning("'groupBy' not provided. Defaulting to 'type'.", call. = FALSE)
      groupBy <- "type"
    }
  }

  if (!groupBy %in% names(group_map)) {
    stop("Invalid 'groupBy'. Choose from: ", paste(names(group_map), collapse=", "))
  }
  col_name <- group_map[[groupBy]]

  # --- 3. Handle Item Selection ---
  if (is.null(selectIds) && interactive()) {
    items_vector <- as.character(interactionData[[col_name]])
    items_vector <- items_vector[!is.na(items_vector)]
    item_counts <- sort(table(items_vector), decreasing = TRUE)
    available_items <- names(item_counts)

    message(if (highlightSelected) "\n--- Select items to HIGHLIGHT ---\n" else
      "\n--- Select items to FILTER ---\n")
    selectIds <- .selectItemsPaged(available_items, groupName = groupBy)
  }

  # --- 4. Prepare Data (Filter/Highlight) ---
  if (highlightSelected) {
    if (is.null(selectIds)) stop("To use 'highlightSelected', you must provide 'selectIds'.")
    pathwayData <- interactionData[interactionData[[col_name]] %in% selectIds, ]
    plotData <- interactionData
  } else {
    if (!is.null(selectIds)) {
      plotData <- interactionData[interactionData[[col_name]] %in% selectIds, ]
    } else {
      plotData <- interactionData
    }
    pathwayData <- NULL
  }

  if (nrow(plotData) == 0) stop("No data to plot after filtering.")

  # --- 5. Prepare Sankey Data ---
  GOprot <- data.frame(
    source = plotData$intersection,
    target = plotData$term_name,
    stringsAsFactors = FALSE
  )
  LNCprot <- data.frame(
    source = plotData$lncRNAId,
    target = plotData$intersection,
    stringsAsFactors = FALSE
  )

  link <- rbind(GOprot, LNCprot)
  link$value <- 1

  node <- data.frame(
    name = unique(c(link$source, link$target)),
    stringsAsFactors = FALSE
  )

  node <- node[!is.na(node$name), , drop = FALSE]

  link$IDsource <- match(link$source, node$name) - 1
  link$IDtarget <- match(link$target, node$name) - 1

  valid_links <- !is.na(link$IDsource) & !is.na(link$IDtarget)
  link <- link[valid_links, ]

  # --- 6. Set Colors ---
  n_nodes <- nrow(node)

  if (is.null(color)) {
    seed_colors <- c("#ff0000", "#00ff00", "#0000ff")
    if (n_nodes > 0) {
      if (n_nodes > length(seed_colors)) {
        mypal <- tryCatch(
          Polychrome::createPalette(n_nodes, seed_colors, M = 1000),
          error = function(e) {
            grDevices::rainbow(n_nodes)
          }
        )
        node$color <- as.character(mypal)
      } else {
        node$color <- seed_colors[1:n_nodes]
      }
    }
  } else {
    node$color <- color
  }

  if (highlightSelected && !is.null(pathwayData)) {
    pathway_nodes <- unique(c(
      pathwayData$intersection,
      pathwayData$term_name,
      pathwayData$lncRNAId
    ))
    node$color[!node$name %in% pathway_nodes] <- "gray"
  }

  link$color <- node$color[link$IDtarget + 1]

  if (highlightSelected && !is.null(pathwayData)) {
    valid_LNC_Int <- paste(pathwayData$lncRNAId, pathwayData$intersection)
    valid_Int_Term <- paste(pathwayData$intersection, pathwayData$term_name)

    current_links_str <- paste(link$source, link$target)
    valid_links_str <- c(valid_LNC_Int, valid_Int_Term)

    link$color[!current_links_str %in% valid_links_str] <- "gray"
  }

  # --- 7. Create Plot ---
  plot_node <- list(
    label = if (showLabels) node$name else NULL,
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
  )

  layout_args <- list(font = list(size = 10))
  if (!is.null(title)) {
    layout_args$title <- title
  }

  fig <- do.call(plotly::layout, c(list(p = fig), layout_args))

  return(fig)
}

# --- Internal Helper Functions  ---

#' @noRd
.selectGroupByInteractive <- function() {
    options <- c("lncRNAId", "target", "term", "type")
    header <- "\n--- Select Grouping Variable ---\nTip: You can skip this menu by providing the 'groupBy' argument directly.\n      e.g., plotInteractions(..., groupBy = \"lncRNAId\")\n"

    menu_lines <- sprintf("%d: \"%s\"", seq_along(options), options)
    menu_text <- paste(menu_lines, collapse = "\n")

    message(header)
    message(menu_text)

    input <- readline("Enter selection number (or 0 to cancel): ")
    input <- trimws(input)

    if (!grepl("^[0-9]+$", input)) {
        message("Invalid selection (not a number).\n")
        return(NULL)
    }

    idx <- as.integer(input)

    if (idx == 0) {
        message("Cancelled.\n")
        return(NULL)
    }

    if (idx < 1 || idx > length(options)) {
        message("Invalid selection (number out of range).\n")
        return(NULL)
    }

    return(options[idx])
}

#' @noRd
.selectItemsPaged <- function(items, groupName, pageSize = 20) {
    n_items <- length(items)
    start_idx <- 1
    selected_items <- character(0)

    while (TRUE) {
        end_idx <- min(start_idx + pageSize - 1, n_items)
        current_page_items <- items[start_idx:end_idx]

        message(sprintf("\n--- Select %s (Showing %d-%d of %d) ---",
                        groupName, start_idx, end_idx, n_items))
        display_indices <- seq_along(current_page_items)
        menu_lines <- sprintf("%d: %s", display_indices, current_page_items)
        message(paste(menu_lines, collapse = "\n"))

        message("\nOptions:")
        message(" - Enter numbers separated by commas to select (e.g., '1,3')")
        message(" - Enter 'a' to select ALL available items")
        if (end_idx < n_items) {
            message(" - Press [Enter] to see the next page")
        }
        message(" - Enter 'q' to finish selection with current choices (if any) or quit")

        input <- readline("Your choice: ")
        input <- trimws(input)

        if (input == "") {
            if (end_idx < n_items) {
                start_idx <- end_idx + 1
                next
            } else {
                message("End of list.\n")
                break
            }
        } else if (tolower(input) == "a") {
            message("Selected ALL items.\n")
            return(items)
        } else if (tolower(input) == "q") {
            break
        } else {
            parts <- strsplit(input, ",")[[1]]
            parts <- trimws(parts)

            valid_parts <- parts[grepl("^[0-9]+$", parts)]

            if (length(valid_parts) < length(parts)) {
                message("Note: Some inputs were not numbers and were ignored.")
            }

            indices <- as.integer(valid_parts)

            valid_indices <- indices[indices >= 1 & indices <= length(current_page_items)]

            if (length(valid_indices) > 0) {
                new_picks <- current_page_items[valid_indices]
                selected_items <- c(selected_items, new_picks)
                message(sprintf("Added %d items. Total selected: %d",
                                length(new_picks), length(unique(selected_items))))
            }

            cont <- readline("Continue selecting? (y/n): ")
            if (tolower(cont) != "y") break
        }
    }

    if (length(selected_items) == 0) {
        message("No specific items selected.\n")
        return(NULL)
    }

    return(unique(selected_items))
}
