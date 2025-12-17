#' Prepare Cis-Interaction Table for Sankey Diagram Plotting
#'
#' Merges gProfiler enrichment results for cis-regulatory interactions
#' with an original table of all cis interactions. This function prepares
#' data for visualization with Sankey diagrams.
#'
#' @param cisGprof data.frame. Results from a gProfiler enrichment analysis
#'   performed on genes from cis-regulatory interactions. Must contain the columns
#'   `term_name`, `intersection`, and `source`.
#' @param cisTable data.frame. The original data frame of cis-regulatory
#'   interactions. Must contain columns `lncRNA_gene` and `partnerRNA_gene`.
#'
#' @return A data.frame suitable for creating Sankey diagrams. The returned
#'   data frame links functional terms to the original cis-interaction partners.
#'   Key columns include:
#'   \itemize{
#'     \item{\code{term_name}}{: Functional term name from gProfiler.}
#'     \item{\code{intersection}}{: Gene ID, corresponding to `targetRNAId`.}
#'     \item{\code{source}}{: Source database of the functional term.}
#'     \item{\code{type}}{: Interaction type, set to "cis".}
#'     \item{\code{lncRNAId}}{: Identifier for the lncRNA gene.}
#'     \item{\code{targetRNAId}}{: Identifier for the target RNA gene.}
#'   }
#'
#' @importFrom dplyr select mutate rename inner_join distinct
#' @importFrom tidyr separate_rows
#' @export
#' @examples
#' # --- 1. Create sample input data frames ---
#'
#' # Sample gProfiler results
#' sampleCisGprof <- data.frame(
#'   term_name = c("GO:0006355", "GO:0006351", "KEGG:04110"),
#'   intersection = c("GENE_A", "GENE_B,GENE_C", "GENE_C"),
#'   source = c("GO:BP", "GO:BP", "KEGG"),
#'   p_value = c(0.01, 0.05, 0.02),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Sample original cis interactions table
#' sampleCisTable <- data.frame(
#'   lncRNA_gene = c("LNC1", "LNC2", "LNC3"),
#'   partnerRNA_gene = c("GENE_A", "GENE_B", "GENE_C"),
#'   distance = c(1000, 2000, 3000),
#'   stringsAsFactors = FALSE
#' )
#'
#' # --- 2. Run the cisInteractions function ---
#' processedCisInteractions <- cisInteractions(
#'   cisGprof = sampleCisGprof,
#'   cisTable = sampleCisTable
#' )
#'
#' print(processedCisInteractions)
#'
#' # The output 'processedCisInteractions' is now ready for plotting functions
#' # like plot_by_terms.
#'
cisInteractions <- function(cisGprof, cisTable) {

  # --- 1. Input Validation ---
  stopifnot(
    "cisGprof must be a data.frame" = is.data.frame(cisGprof),
    "cisTable must be a data.frame" = is.data.frame(cisTable)
  )
  requiredGprofCols <- c("term_name", "intersection", "source")
  requiredCisCols <- c("lncRNA_gene", "partnerRNA_gene")
  stopifnot(
    "cisGprof is missing required columns" = all(requiredGprofCols %in% names(cisGprof)),
    "cisTable is missing required columns" = all(requiredCisCols %in% names(cisTable))
  )

  # --- 2. Process Data using dplyr/tidyr pipeline ---
  processedGprof <- cisGprof |>
    dplyr::select(term_name, intersection, source) |>
    dplyr::mutate(type = "cis") |>
    tidyr::separate_rows(intersection, sep = ",")

  processedCisTable <- cisTable |>
    dplyr::rename(
      lncRNAId = lncRNA_gene,
      targetRNAId = partnerRNA_gene
    )

  # Merge the two data frames
  finalTable <- dplyr::inner_join(
    processedGprof,
    processedCisTable,
    by = c("intersection" = "targetRNAId")
  ) |>
    dplyr::rename(targetRNAId = intersection) |>
    dplyr::distinct()

  return(as.data.frame(finalTable))
}
