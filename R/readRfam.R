#' Read and Filter an Rfam hmmscan Result File
#'
#' This function reads a `domtblout` file from an hmmscan search against the
#' Rfam database, filters the results by an e-value cutoff, and returns a unique
#' list of query names that pass the filter.
#'
#' @param rfamOutfile A character string specifying the path to the hmmscan
#'   `domtblout` output file.
#' @param evalCutoff A numeric value for the sequence e-value cutoff. Hits with
#'   an e-value less than this threshold will be retained. Defaults to 1e-2.
#'
#' @return A unique character vector of query names that pass the e-value
#'   filter.
#'
#' @keywords rfam lncRNA hmmscan rhmmer
#' @importFrom rhmmer read_domtblout
#' @export
#' @examples
#' # --- 1. Create a temporary dummy domtblout file for a reproducible example ---
#' tempRfamFile <- tempfile()
#'
#' # Create content mimicking the space-delimited domtblout format.
#' # We only need to populate columns 4 (query_name) and 12 (sequence_evalue).
#' filler <- paste(rep("-", 20), collapse = " ")
#' fileContent <- c(
#'   "# Dummy Rfam domtblout file for example",
#'   paste("RF00005", "-", "-", "Seq_A", filler, "1.0e-5", "-"), # Should pass
#'   paste("RF00017", "-", "-", "Seq_A", filler, "2.0e-4", "-"), # Should pass (duplicate query)
#'   paste("RF00001", "-", "-", "Seq_B", filler, "9.0e-1", "-"), # Should be filtered out
#'   paste("RF00002", "-", "-", "Seq_C", filler, "5.0e-3", "-")  # Should pass
#' )
#'
#' # Write the content to the temporary file
#' writeLines(fileContent, tempRfamFile)
#'
#' # --- 2. Run the readRfam function ---
#' if (requireNamespace("rhmmer", quietly = TRUE)) {
#'   passingQueries <- readRfam(rfamOutfile = tempRfamFile, evalCutoff = 1e-2)
#'   print("Unique query names passing the e-value cutoff:")
#'   print(passingQueries) # Expected output: c("Seq_A", "Seq_C")
#' }
#'
#' # --- 3. Clean up the temporary file ---
#' unlink(tempRfamFile)
#'
readRfam <- function(rfamOutfile, evalCutoff = 1e-2) {

  # --- 1. Input Validation ---
  stopifnot(
    "'rfamOutfile' must be a single character string" =
      (is.character(rfamOutfile) && length(rfamOutfile) == 1),
    "File specified in 'rfamOutfile' does not exist" =
      file.exists(rfamOutfile)
  )

  # --- 2. Read and Process Data ---
  rfamHits <- rhmmer::read_domtblout(file = rfamOutfile)

  if (nrow(rfamHits) == 0) {
    return(character(0))
  }

  # --- 3. Filter, Extract, and Unify ---
  # This restores the intended filtering logic
  passingQueries <- rfamHits[rfamHits$sequence_evalue < evalCutoff, "query_name"]
  uniqueNames <- unique(passingQueries)

  return(uniqueNames)
}
