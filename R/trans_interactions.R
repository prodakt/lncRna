#' Prepare Trans-Interaction Table for Sankey Diagram Plotting
#'
#' Processes and merges gProfiler enrichment results for trans-regulatory interactions with an original
#' table of all trans interactions. This function mirrors `cis_interactions` but is tailored for
#' trans-interactions, preparing data for Sankey diagram visualization, especially using the
#' `plot_by_terms` function.
#'
#' @param trans_gprof data.frame. A data frame resulting from a gProfiler enrichment analysis on genes
#'        involved in trans-regulatory interactions. This should be the direct output from gProfiler,
#'        containing columns like `term_name` (functional term name), `intersection` (genes or regions),
#'        and `source` (term database source).
#' @param trans_table data.frame. The original data frame detailing all identified trans-regulatory interactions.
#'        This table is expected to have columns for `lncRNA.id` (lncRNA gene identifier) and
#'        `targetRNA.id` (protein-coding or target RNA gene identifier) for trans-interactions.
#'
#' @return data.frame.
#' Returns a processed and merged data frame designed for creating Sankey diagrams of trans-regulatory
#' interactions. It combines functional term enrichment from `trans_gprof` with trans-interaction details
#' from `trans_table`.  The returned data frame links enriched functional terms (`term_name`, `source`, `type`)
#' to their corresponding trans-interaction partners (`lncRNA.id`, `targetRNA.id`) through the `intersection`
#' column (representing `targetRNA.id` from `trans_table`).
#' The output data frame includes these essential columns:
#'   \itemize{
#'     \item{\code{term_name}}{: Functional term name from gProfiler.}
#'     \item{\code{intersection}}{: Gene or genomic region intersection, corresponding to `targetRNA.id` in `trans_table`.}
#'     \item{\code{source}}{: Source database of the functional term (from gProfiler).}
#'     \item{\code{type}}{: Interaction type, set to "trans" for all entries.}
#'     \item{\code{lncRNA.id}}{: Identifier for the lncRNA gene in the trans-interaction (from `trans_table`).}
#'     \item{\code{targetRNA.id}}{: Identifier for the protein-coding or target RNA gene (from `trans_table`).}
#'   }
#'
#' @examples
#' # Assuming 'trans_gProfiler_results' and 'original_trans_interactions' are your input data frames
#' # (replace with your actual data frame names)
#'
#' # Process trans-interactions for Sankey diagram plotting
#' # processed_trans_interactions <- trans_interactions(trans_gprof = trans_gProfiler_results,
#' #                                                  trans_table = original_trans_interactions)
#'
#' # Use the processed table with plot_by_terms to generate a Sankey diagram:
#' # sankey_plot_trans <- plot_by_terms(data = processed_trans_interactions,
#' #                                     title = "Sankey Diagram of Trans-Regulatory Interactions")
#' # sankey_plot_trans # Display the plot
#'
#' @import dplyr
#' @import tidyr
#' @export
trans_interactions <- function(trans_gprof, trans_table) {
  trans_gprof$type <- "trans"
  trans_gprof <- trans_gprof[,c("term_name", "intersection","source","type")]
  trans_gprof <- separate_rows(trans_gprof, "intersection", sep = ",")
  trans_gprof <- as.data.frame(trans_gprof)
  trans_table_final <- merge(trans_gprof, trans_table[,c("lncRNA.id", "targetRNA.id")], by.x = "intersection",
                             by.y = "targetRNA.id", all.x = T)
  trans_table_final <- unique(trans_table_final)
  return(trans_table_final)
}
#######################################
# Helper Functions Documentation
