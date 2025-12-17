#' Prepare Data Sets for Performance Evaluation
#'
#' Filters coding potential results to include only sequences present in
#' provided test sets. It annotates each sequence as non-coding (nc) or
#' protein-coding (cds), creating a summary object ready for various evaluation
#' functions (e.g., `BestTool`, `BestToolAtleast`).
#'
#' @param codPotList A list object, typically from `aggregateCodPot()`. Must
#'   contain `$seqIDs` (character vector) and `$tools` (a list of
#'   numeric vectors).
#' @param ncTest A character vector of sequence IDs known to be non-coding.
#' @param cdsTest A character vector of sequence IDs known to be protein-coding.
#'
#' @return A list containing elements filtered to include only sequences from
#'   `ncTest` or `cdsTest`:
#'   \item{seqIDs}{Filtered character vector of sequence IDs.}
#'   \item{tools}{Filtered list of tool prediction vectors.}
#'   \item{type}{Character vector with annotation ("nc" or "cds").}
#'   \item{isNC}{Integer vector: 1 if type is "nc", 0 if "cds".}
#'   \item{sums}{Integer vector with the sum of predictions across all tools.}
#' @export
#'
#' @examples
#' # --- 1. Create mock data mimicking package outputs ---
#' # a) Output from aggregateCodPot()
#' mockCodPotList <- list(
#'   seqIDs = c("nc_seq1", "cds_seq1", "nc_seq2", "other_seq"),
#'   tools = list(
#'     CPC2 = c(1, 0, 1, 1),
#'     PLEK = c(1, 1, 0, 0)
#'   )
#' )
#' 
#' # b) Outputs from createTrainTestSets()
#' mockNcSets <- list(
#'   nc.train = c("nc_seq3", "nc_seq4"),
#'   nc.test = c("nc_seq1", "nc_seq2")
#' )
#' mockCdsSets <- list(
#'   cds.train = c("cds_seq2"),
#'   cds.test = c("cds_seq1")
#' )
#'
#' # --- 2. Run the function to prepare the evaluation set ---
#' evaluationSummary <- prepareEvaluationSets(
#'   codPotList = mockCodPotList,
#'   ncTest = mockNcSets$nc.test,
#'   cdsTest = mockCdsSets$cds.test
#' )
#'
#' # --- 3. Inspect the prepared data ---
#' # Note: "other_seq" was filtered out as it was not in the test sets.
#' print(evaluationSummary)
#'
prepareEvaluationSets <- function(codPotList, ncTest, cdsTest) {
    # --- Input Validation ---
    stopifnot(
        "'codPotList' must be a list with 'seqIDs' and 'tools' elements" =
            is.list(codPotList) && all(c("seqIDs", "tools") %in% names(codPotList)),
        "'codPotList$tools' must be a list" = is.list(codPotList$tools),
        "'ncTest' must be a character vector" = is.character(ncTest),
        "'cdsTest' must be a character vector" = is.character(cdsTest)
    )
    
    original_seqIDs <- codPotList$seqIDs
    original_tools <- codPotList$tools
    
    if (length(original_seqIDs) == 0) {
        warning("Input 'codPotList' contains no sequence IDs.")
        return(list(seqIDs = character(0), tools = list(), type = character(0),
                    isNC = integer(0), sums = integer(0)))
    }
    
    # --- Identify Indices to Keep ---
    indices_to_keep <- which(original_seqIDs %in% c(ncTest, cdsTest))
    
    if (length(indices_to_keep) == 0) {
        warning("None of the sequences in codPotList were found in the test sets.")
        return(list(seqIDs = character(0), tools = list(), type = character(0),
                    isNC = integer(0), sums = integer(0)))
    }
    
    # --- Filter Data ---
    filtered_seqIDs <- original_seqIDs[indices_to_keep]
    filtered_tools <- lapply(original_tools, function(tool_vector) {
        if (length(tool_vector) != length(original_seqIDs)) {
            warning("Length mismatch for a tool vector. Skipping this tool.")
            return(NULL)
        }
        tool_vector[indices_to_keep]
    })
    filtered_tools <- Filter(Negate(is.null), filtered_tools)
    
    # --- Create Annotations for Filtered Data ---
    n_filtered <- length(filtered_seqIDs)
    filtered_type <- character(n_filtered)
    filtered_type[filtered_seqIDs %in% ncTest] <- "nc"
    filtered_type[filtered_seqIDs %in% cdsTest] <- "cds"
    
    filtered_isNC <- as.integer(filtered_type == "nc")
    
    # --- Calculate Sums for Filtered Data ---
    if (length(filtered_tools) > 0) {
        filtered_tool_matrix <- do.call(cbind, filtered_tools)
        filtered_sums <- rowSums(filtered_tool_matrix, na.rm = TRUE)
    } else {
        warning("No valid tool data found. 'sums' will be zero.")
        filtered_sums <- integer(n_filtered)
    }
    
    # --- Assemble the Final List ---
    result_list <- list(
        seqIDs = filtered_seqIDs,
        tools = filtered_tools,
        type = filtered_type,
        isNC = filtered_isNC,
        sums = as.integer(filtered_sums)
    )
    
    return(result_list)
}