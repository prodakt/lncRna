#' Read and Filter a Pfam hmmscan Result File
#'
#' This function reads a `domtblout` file from an hmmscan search against the
#' Pfam database, filters the results by an e-value cutoff, and returns a unique
#' list of query names.
#'
#' @param pfamOutfile A character string specifying the path to the hmmscan
#'   `domtblout` output file.
#' @param evalCutoff A numeric value for the sequence e-value cutoff. Hits with
#'   an e-value less than this threshold will be retained. Defaults to 1e-2.
#'
#' @return A unique character vector of query names (with protein-specific
#'   suffixes like '.p1' removed) that pass the e-value filter.
#'
#' @keywords pfam lncRNA hmmscan rhmmer
#' @importFrom rhmmer read_domtblout
#' @export
#' @examples
#' # --- 1. Create a temporary dummy domtblout file for a reproducible example ---
#' tempPfamFile <- tempfile()
#'
#' # Create content mimicking the space-delimited domtblout format.
#' # We only need to populate columns 4 (query_name) and 12 (sequence_evalue).
#' # The other 20 columns can be placeholders.
#' filler <- paste(rep("-", 20), collapse = " ")
#' fileContent <- c(
#'   "# Dummy domtblout file for example",
#'   paste("target1", "-", "-", "Seq1.p1", filler, "1.5e-5", "-"), # Should pass
#'   paste("target2", "-", "-", "Seq1.p2", filler, "2.0e-4", "-"), # Should pass (duplicate after cleaning)
#'   paste("target3", "-", "-", "Seq2.p1", filler, "5.0e-2", "-"), # Should be filtered out
#'   paste("target4", "-", "-", "Seq3",    filler, "1.0e-10", "-") # Should pass (no suffix)
#' )
#'
#' # Write the content to the temporary file
#' writeLines(fileContent, tempPfamFile)
#'
#' # --- 2. Run the readPfam function ---
#' if (requireNamespace("rhmmer", quietly = TRUE)) {
#'   passingQueries <- readPfam(pfamOutfile = tempPfamFile, evalCutoff = 1e-2)
#'   print("Unique query names passing the e-value cutoff:")
#'   print(passingQueries) # Expected output: c("Seq1", "Seq3")
#' }
#'
#' # --- 3. Clean up the temporary file ---
#' unlink(tempPfamFile)
#'
readPfam <- function(pfamOutfile, evalCutoff = 1e-2) {

  # --- 1. Input Validation ---
  stopifnot(
    "'pfamOutfile' must be a single character string" =
      (is.character(pfamOutfile) && length(pfamOutfile) == 1),
    "File specified in 'pfamOutfile' does not exist" =
      file.exists(pfamOutfile)
  )

  # --- 2. Read and Process Data ---
  # read_domtblout from rhmmer handles the complex file format
  pfamHits <- rhmmer::read_domtblout(file = pfamOutfile)

  # Return empty vector if no hits were found in the file
  if (nrow(pfamHits) == 0) {
    return(character(0))
  }

  # Filter by e-value, extract query names, clean, and find unique
  passingQueries <- pfamHits[pfamHits$sequence_evalue < evalCutoff, "query_name"]
  cleanedNames <- gsub("\\.p.*$", "", passingQueries)
  uniqueNames <- unique(cleanedNames)

  return(uniqueNames)
}
