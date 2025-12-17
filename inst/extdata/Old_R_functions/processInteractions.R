#' Process and Merge Genomic Interactions
#'
#' Processes genomic profiler data and merges it with an interaction table to
#' create a unified dataset. Allows for custom column names for lncRNA and
#' target identifiers.
#'
#' @param gprof A data.frame from gProfiler. Must include `term_name` and
#'   `intersection` columns.
#' @param interactionTable A data.frame with interaction data. Must include
#'   columns for lncRNA and target identifiers.
#' @param type A character string specifying the interaction type (e.g., "cis").
#' @param lncRnaCol An optional character string specifying the column name in
#'   `interactionTable` for lncRNA IDs. Used if `lncRNAId` is not present.
#' @param targetCol An optional character string specifying the column name in
#'   `interactionTable` for target IDs. Used if `targetRNAId` is not present.
#'
#' @return A unique data.frame merging the processed gProfiler data and
#'   interaction data.
#'
#' @importFrom tidyr separate_rows
#' @importFrom dplyr select rename left_join distinct any_of
#' @export
#' @examples
#' # Example 1: Using default column names
#' gprofData <- data.frame(
#'   term_name = "response to stress",
#'   intersection = "GENE_A",
#'   source = "GO:BP"
#' )
#' intTable <- data.frame(lncRNAId = "LNC_1", targetRNAId = "GENE_A")
#' result1 <- processInteractions(gprofData, intTable, type = "trans")
#' print(result1)
#'
#' # Example 2: Using custom column names
#' gprofData2 <- data.frame(
#'   term_name = "cell activation",
#'   intersection = "GENE_B",
#'   source = "GO:BP"
#' )
#' intTable2 <- data.frame(Query = "LNC_2", Target = "GENE_B")
#' result2 <- processInteractions(
#'   gprofData2,
#'   intTable2,
#'   type = "trans",
#'   lncRnaCol = "Query",
#'   targetCol = "Target"
#' )
#' print(result2)
processInteractions <- function(gprof, interactionTable, type,
                                lncRnaCol = NULL, targetCol = NULL) {
  # --- 1. Validate Inputs ---
  stopifnot(
    "'gprof' must be a data.frame" = is.data.frame(gprof),
    "'interactionTable' must be a data.frame" = is.data.frame(interactionTable),
    "'gprof' is missing required columns" = all(c("term_name", "intersection") %in% names(gprof)),
    "'type' must be a single, non-empty string" = (is.character(type) && length(type) == 1 && nzchar(type))
  )

  # --- 2. Standardize Column Names in interactionTable ---
  if ("lncRNAId" %in% names(interactionTable) && "targetRNAId" %in% names(interactionTable)) {
    if (!is.null(lncRnaCol) || !is.null(targetCol)) {
      warning("Default columns 'lncRNAId' and 'targetRNAId' found; ignoring 'lncRnaCol' and 'targetCol'.")
    }
  } else {
    stopifnot(
      "'lncRnaCol' and 'targetCol' must be specified when default columns are missing" =
        (!is.null(lncRnaCol) && !is.null(targetCol)),
      "Specified 'lncRnaCol' not found in interactionTable" = (lncRnaCol %in% names(interactionTable)),
      "Specified 'targetCol' not found in interactionTable" = (targetCol %in% names(interactionTable))
    )
    interactionTable <- interactionTable |>
      dplyr::rename(lncRNAId = !!lncRnaCol, targetRNAId = !!targetCol)
  }

  # --- 3. Process and Merge Data using a dplyr/tidyr Pipeline ---
  finalTable <- gprof |>
    dplyr::select(dplyr::any_of(c("term_name", "intersection", "source"))) |>
    dplyr::mutate(type = type) |>
    tidyr::separate_rows(intersection, sep = ",") |>
    dplyr::left_join(
      interactionTable |> dplyr::select(lncRNAId, targetRNAId),
      by = c("intersection" = "targetRNAId")
    ) |>
    dplyr::distinct()

  return(as.data.frame(finalTable))
}
