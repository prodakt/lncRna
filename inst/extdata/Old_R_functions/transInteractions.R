#' Prepare Trans-Interaction Table for Sankey Diagram Plotting
#'
#' Merges gProfiler enrichment results for trans-regulatory interactions with an
#' original table of all trans interactions. This function prepares data for
#' visualization with Sankey diagrams.
#'
#' @param transGprof A data.frame from a gProfiler enrichment analysis on genes
#'   from trans-regulatory interactions. Must contain `term_name`,
#'   `intersection`, and `source` columns.
#' @param transTable A data.frame of trans-regulatory interactions. Must
#'   contain `lncRNAId` and `targetRNAId` columns.
#'
#' @return A data.frame suitable for creating Sankey diagrams, linking
#'   functional terms to the original trans-interaction partners. Key columns
#'   include `term_name`, `targetRNAId`, `source`, `type`, and `lncRNAId`.
#'
#' @importFrom dplyr select mutate left_join distinct rename
#' @importFrom tidyr separate_rows
#' @export
#' @examples
#' # --- 1. Create sample input data frames ---
#'
#' # Sample gProfiler results for trans-interactions
#' sampleTransGprof <- data.frame(
#'   term_name = c("cell cycle", "DNA repair"),
#'   intersection = c("GENE_X,GENE_Y", "GENE_Z"),
#'   source = c("KEGG", "GO:BP"),
#'   p_value = c(0.01, 0.03),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Sample original trans-interactions table
#' sampleTransTable <- data.frame(
#'   lncRNAId = c("LNC_A", "LNC_B", "LNC_C"),
#'   targetRNAId = c("GENE_X", "GENE_Y", "GENE_Z"),
#'   rValue = c(0.8, -0.9, 0.85),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Run the transInteractions function ---
#' processedTransInteractions <- transInteractions(
#'   transGprof = sampleTransGprof,
#'   transTable = sampleTransTable
#' )
#'
#' print(processedTransInteractions)
#'
#' # The output 'processedTransInteractions' is now ready for plotting functions.
#'
transInteractions <- function(transGprof, transTable) {

  # --- 1. Input Validation ---
  stopifnot(
    "'transGprof' must be a data.frame" = is.data.frame(transGprof),
    "'transTable' must be a data.frame" = is.data.frame(transTable)
  )
  requiredGprofCols <- c("term_name", "intersection", "source")
  requiredTransCols <- c("lncRNAId", "targetRNAId")
  stopifnot(
    "'transGprof' is missing required columns" = all(requiredGprofCols %in% names(transGprof)),
    "'transTable' is missing required columns" = all(requiredTransCols %in% names(transTable))
  )

  # --- 2. Process Data using a dplyr/tidyr pipeline ---
  finalTable <- transGprof |>
    dplyr::select(term_name, intersection, source) |>
    dplyr::mutate(type = "trans") |>
    tidyr::separate_rows(intersection, sep = ",") |>
    dplyr::left_join(
      transTable |> dplyr::select(lncRNAId, targetRNAId),
      by = c("intersection" = "targetRNAId")
    ) |>
    dplyr::rename(targetRNAId = intersection) |>
    dplyr::distinct()

  return(as.data.frame(finalTable))
}
