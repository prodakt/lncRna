#' Prepare Cis-Interaction Table for Sankey Diagram Plotting
#'
#' Creates a processed data table by merging gProfiler enrichment results for cis-regulatory interactions
#' with an original table of all cis interactions. This function is specifically designed to prepare
#' data for visualization with Sankey diagrams, especially using the `plot_by_terms` function.
#'
#' @param cis_gprof data.frame. A data frame containing the results of a gProfiler enrichment analysis
#'        performed on genes involved in cis-regulatory interactions. This table should be the direct output
#'        from running gProfiler on a list of genes identified in cis-interactions. It must contain columns
#'        such as `term_name` (functional term name), `intersection` (genes or genomic regions associated with the term),
#'        and `source` (the data source of the term).
#' @param cis_table data.frame. The original data frame detailing all identified cis-regulatory interactions.
#'        This table should contain at least columns for `lncRNA_gene` (identifier for the lncRNA gene)
#'        and `partnerRNA_gene` (identifier for the protein-coding gene or target RNA gene) involved in the cis-interaction.
#'
#' @return data.frame.
#' Returns a merged and processed data frame suitable for creating Sankey diagrams of cis-regulatory interactions.
#' The returned data frame includes columns from both `cis_gprof` and `cis_table`, effectively linking
#' functional terms enriched in gProfiler (`term_name`, `source`, `type`) to the original cis-interaction
#' partners (`lncRNA.id`, `targetRNA.id`) via the `intersection` (which represents the `targetRNA.id` from `cis_table`).
#' The output data frame contains the following key columns:
#'   \itemize{
#'     \item{\code{term_name}}{: Functional term name from gProfiler.}
#'     \item{\code{intersection}}{: Gene or genomic region intersection, corresponding to `targetRNA.id` in `cis_table`.}
#'     \item{\code{source}}{: Source database of the functional term (from gProfiler).}
#'     \item{\code{type}}{: Interaction type, set to "cis" for all entries in this table.}
#'     \item{\code{lncRNA.id}}{: Identifier for the lncRNA gene involved in the cis-interaction (from `cis_table`).}
#'     \item{\code{targetRNA.id}}{: Identifier for the protein-coding or target RNA gene (from `cis_table`).}
#'   }
#'
#' @examples
#' # Assuming 'cis_gProfiler_results' and 'original_cis_interactions' are your input data frames
#' # (replace with your actual data frame names)
#'
#' # Prepare the cis-interactions table for Sankey diagram plotting
#' # processed_cis_interactions <- cis_interactions(cis_gprof = cis_gostres,
#' #                                                cis_table = cis_table)
#'
#' # Now 'processed_cis_interactions' can be used as input for plot_by_terms function:
#' # sankey_plot_cis <- plot_by_terms(data = processed_cis_interactions,
#' #                                   title = "Sankey Diagram of Cis-Regulatory Interactions")
#' # sankey_plot_cis # Display the plot
#'
#' @import dplyr
#' @import tidyr
#' @export
cis_interactions <- function(cis_gprof, cis_table) {
  cis_gprof$type <- "cis"
  cis_gprof <-    cis_gprof[,c("term_name", "intersection","source","type")]
  cis_gprof <- separate_rows(cis_gprof, "intersection", sep = ",")
  cis_gprof <- as.data.frame(cis_gprof)
  names(cis_table)[names(cis_table) == 'lncRNA_gene'] <- 'lncRNA.id'
  names(cis_table)[names(cis_table) == 'partnerRNA_gene'] <- 'targetRNA.id'
  cis_table_final <- merge(cis_gprof, cis_table[,c("lncRNA.id", "targetRNA.id")], by.x = "intersection",
                           by.y = "targetRNA.id", all.x = T)
  cis_table_final <- unique(cis_table_final)
  return(cis_table_final)
}
