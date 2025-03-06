#' Process and Merge Genomic Interactions
#'
#' Processes genomic profiler data and merges it with an interaction table to create a unified
#' dataset of interactions. Requires specification of interaction type and ensures correct column
#' naming for lncRNA and target identifiers.
#'
#' @param gprof data.frame. Genomic profiler data (e.g., output from gProfiler).
#'        Must include at least the columns "term_name" and "intersection", and optionally "source".
#' @param interaction_table data.frame. Interaction data with columns for lncRNA and target identifiers.
#'        Must include either "lncRNA.id" and "targetRNA.id", or columns specified by `lncRNA_col` and `target_col`.
#' @param type character. The type of interaction to assign (e.g., "cis" or "trans"). Must be specified.
#' @param lncRNA_col character, optional. Name of the column in `interaction_table` containing lncRNA identifiers.
#'        Required if "lncRNA.id" is not present in `interaction_table`. Defaults to NULL.
#' @param target_col character, optional. Name of the column in `interaction_table` containing target identifiers.
#'        Required if "targetRNA.id" is not present in `interaction_table`. Defaults to NULL.
#'
#' @return data.frame. A merged data frame containing the processed genomic profiler data and interaction data,
#'         with unique rows. Includes columns such as "intersection", "term_name", "source", "type", and "lncRNA.id".
#'
#' @examples
#' # Example with default column names
#' gprof <- data.frame(term_name = "response to stress", intersection = "ENSMUSG00000000093", source = "GO:BP")
#' int_table <- data.frame(lncRNA.id = "ENSMUSG00000106858", targetRNA.id = "ENSMUSG00000000093")
#' result <- process_interactions(gprof, int_table, type = "trans")
#'
#' # Example with custom column names
#' gprof2 <- data.frame(term_name = "cell activation", intersection = "ENSMUSG00000000731", source = "GO:BP")
#' int_table2 <- data.frame(Query = c("ENSMUSG00000002769", "ENSMUSG00000089736"), Target = "ENSMUSG00000000731")
#' result2 <- process_interactions(gprof2, int_table2, type = "trans", lncRNA_col = "Query", target_col = "Target")
#'
#' @import tidyr
#' @import dplyr
#' @export
process_interactions <- function(gprof, interaction_table, type, lncRNA_col = NULL, target_col = NULL) {
  # --- Step 1: Validate Inputs ---
  
  # Validate gprof
  if (!is.data.frame(gprof)) {
    stop("'gprof' must be a data frame.")
  }
  required_gprof_cols <- c("term_name", "intersection")
  if (!all(required_gprof_cols %in% colnames(gprof))) {
    stop("'gprof' must contain columns: ", paste(required_gprof_cols, collapse = ", "))
  }
  
  # Validate interaction_table
  if (!is.data.frame(interaction_table)) {
    stop("'interaction_table' must be a data frame.")
  }
  
  # Validate type (mandatory)
  if (missing(type) || is.null(type) || !is.character(type) || nchar(type) == 0) {
    stop("'type' must be a non-empty character string specifying the interaction type (e.g., 'cis' or 'trans' etc....).")
  }
  
  # --- Step 2: Validate and Rename Columns in interaction_table ---
  # Check for default columns "lncRNA.id" and "targetRNA.id"
  has_default_lncRNA <- "lncRNA.id" %in% colnames(interaction_table)
  has_default_target <- "targetRNA.id" %in% colnames(interaction_table)
  
  if (has_default_lncRNA && has_default_target) {
    # Default columns are present; no renaming needed
    if (!is.null(lncRNA_col) || !is.null(target_col)) {
      warning("Default columns 'lncRNA.id' and 'targetRNA.id' are present; ignoring 'lncRNA_col' and 'target_col'.")
    }
  } else {
    # If either default column is missing, require lncRNA_col and target_col
    if (is.null(lncRNA_col) || is.null(target_col)) {
      warning("Missing default columns 'lncRNA.id' and/or 'targetRNA.id' in interaction_table. ",
              "You must specify 'lncRNA_col' and 'target_col'.")
      stop("'lncRNA_col' and 'target_col' must be specified when default columns are missing.")
    }
    
    # Validate specified columns exist in interaction_table
    if (!(lncRNA_col %in% colnames(interaction_table))) {
      stop("lncRNA_col '", lncRNA_col, "' not found in interaction_table.")
    }
    if (!(target_col %in% colnames(interaction_table))) {
      stop("target_col '", target_col, "' not found in interaction_table.")
    }
    
    # Rename specified columns to standard names
    colnames(interaction_table)[colnames(interaction_table) == lncRNA_col] <- "lncRNA.id"
    colnames(interaction_table)[colnames(interaction_table) == target_col] <- "targetRNA.id"
  }
  
  # --- Step 3: Prepare gprof Data ---
  # Add mandatory type column
  gprof$type <- type
  
  # Select relevant columns
  gprof_cols <- c("term_name", "intersection", "source", "type")
  gprof_cols <- gprof_cols[gprof_cols %in% colnames(gprof)]
  gprof <- gprof[, gprof_cols, drop = FALSE]
  
  # Separate intersection column into rows
  gprof <- tidyr::separate_rows(gprof, intersection, sep = ",")
  gprof <- as.data.frame(gprof)
  
  # --- Step 4: Merge Data Frames ---
  table_final <- merge(gprof, interaction_table[, c("lncRNA.id", "targetRNA.id")], 
                       by.x = "intersection", by.y = "targetRNA.id", all.x = TRUE)
  
  # Remove duplicates
  table_final <- unique(table_final)
  
  return(table_final)
}
