#' Annotate Genomic Interactions with Functional Enrichment Results
#'
#' Merges a table of genomic interactions (e.g., lncRNA-mRNA) with functional
#' enrichment results (e.g., from g:Profiler), annotating each interaction
#' with relevant functional terms.
#'
#' @param gostResult A `data.frame` from a g:Profiler analysis (`gost` function),
#'   containing at least "term_name" and "intersection" columns.
#' @param interactionTable A `data.frame` of interactions. Must contain columns
#'   for lncRNA and target identifiers.
#' @param type A character string specifying the interaction type to be added
#'   to the output (e.g., "cis", "trans", "LncTar").
#' @param lncRnaCol An optional character string specifying the column name in
#'   `interactionTable` for lncRNA IDs. Defaults to `lncRNAId`.
#' @param targetCol An optional character string specifying the column name in
#'   `interactionTable` for target RNA IDs. Defaults to `targetRNAId`.
#'
#' @return A `data.frame` where enrichment terms are merged with the interactions
#'   they are associated with.
#' @export
#'
#' @examples
#' # --- 1. Create mock input data ---
#' # a) Mock g:Profiler result
#' mockGostResult <- data.frame(
#'   term_name = c("response to stress", "cell activation"),
#'   intersection = c("TARGET_G1,TARGET_G2", "TARGET_G3"),
#'   source = c("GO:BP", "GO:BP"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # b) Mock interaction table with standard column names
#' mockInteractionTable1 <- data.frame(
#'   lncRNAId = c("LNC_G1", "LNC_G2", "LNC_G3"),
#'   targetRNAId = c("TARGET_G1", "TARGET_G2", "TARGET_G3")
#' )
#'
#' # --- 2. Run with standard column names ---
#' annotated_trans <- annotateInteractions(
#'   gostResult = mockGostResult,
#'   interactionTable = mockInteractionTable1,
#'   type = "trans"
#' )
#' print(annotated_trans)
#'
#' # --- 3. Run with custom column names ---
#' # a) Mock interaction table with custom column names
#' mockInteractionTable2 <- data.frame(
#'   Query = "LNC_G1",
#'   Target = "TARGET_G1"
#' )
#'
#' annotated_lnctar <- annotateInteractions(
#'   gostResult = mockGostResult,
#'   interactionTable = mockInteractionTable2,
#'   type = "LncTar",
#'   lncRnaCol = "Query",
#'   targetCol = "Target"
#' )
#' print(annotated_lnctar)
#'
annotateInteractions <- function(gostResult, interactionTable, type,
                                 lncRnaCol = "lncRNAId",
                                 targetCol = "targetRNAId") {
    
    # --- 1. Input Validation ---
    stopifnot(
        "'gostResult' must be a data.frame" = is.data.frame(gostResult),
        "'interactionTable' must be a data.frame" = is.data.frame(interactionTable),
        "'type' must be a single character string" = is.character(type) && length(type) == 1
    )
    
    required_gost_cols <- c("term_name", "intersection")
    if (!all(required_gost_cols %in% colnames(gostResult))) {
        stop("'gostResult' must contain columns: ",paste(required_gost_cols, collapse = ", "))
    }
    
    if (!(lncRnaCol %in% colnames(interactionTable))) {
        stop("lncRNA column '", lncRnaCol, "' not found in interactionTable.")
    }
    if (!(targetCol %in% colnames(interactionTable))) {
        stop("Target column '", targetCol, "' not found in interactionTable.")
    }
    
    # --- 2. Prepare and expand gostResult  ---
    gostResult$type <- type
    
    # Select and reorder columns, keeping 'source' if it exists
    gost_cols <- c("term_name", "intersection", "source", "type")
    gost_cols_present <- intersect(gost_cols, colnames(gostResult))
    gost_prepared <- gostResult[, gost_cols_present, drop = FALSE]
    
    # Split intersection strings into a list of character vectors
    intersection_list <- strsplit(gost_prepared$intersection, split = ",")
    
    # Calculate how many times each row of gost_prepared needs to be repeated
    times_to_repeat <- vapply(intersection_list, length, integer(1))
    
    # Expand the data frame by repeating rows
    expanded_gost <- gost_prepared[rep(seq_len(nrow(gost_prepared)), times = times_to_repeat), ]
    
    # Replace the comma-separated strings with the unlisted individual gene IDs
    expanded_gost$intersection <- unlist(intersection_list)
    rownames(expanded_gost) <- NULL 
    
    # --- 3. Prepare interactionTable ---
    interactionTable_std <- interactionTable[, c(lncRnaCol, targetCol), drop = FALSE]
    colnames(interactionTable_std) <- c("lncRNAId", "targetRNAId")
    
    # --- 4. Merge Data Frames ---
    final_table <- merge(
        expanded_gost,
        interactionTable_std,
        by.x = "intersection",
        by.y = "targetRNAId",
        all.x = TRUE 
    )
    
    return(unique(final_table))
}