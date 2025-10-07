#' Analyze cis-regulatory interactions from FEELnc output
#'
#' This function reads and filters the 'classes' output file from FEELnc to
#' identify cis-regulatory relationships between lncRNAs and mRNAs based on
#' user-defined criteria.
#'
#' @param feelncClasses character. The path to the "classes" output file from
#'   FEELnc. The file must exist.
#' @param lncRNAs character (optional). A vector of lncRNA gene or transcript IDs
#'   to filter the results. If `NULL` (default), all lncRNAs are considered.
#' @param mRNAs character (optional). A vector of mRNA gene or transcript IDs
#'   to filter the results. If `NULL` (default), all mRNAs are considered.
#' @param isBest logical. If `TRUE` (default), only interactions flagged as the
#'   'best' (isBest == 1) in the FEELnc output are retained.
#' @param lncRnaLevel character. Specifies whether to filter `lncRNAs` at the
#'   "gene" or "transcript" level. Defaults to "transcript".
#' @param mRnaLevel character. Specifies whether to filter `mRNAs` at the
#'   "gene" or "transcript" level. Defaults to "gene".
#' @param maxDist numeric. The maximum allowed distance in base pairs for an
#'   interaction to be considered 'cis'. Defaults to 100000.
#'
#' @return A data.frame containing the filtered cis-interactions. The columns
#'   'lncRNA_gene' and 'partnerRNA_gene' are renamed to 'lncRNAId' and
#'   'targetRNAId' respectively for clarity.
#'
#' @keywords cis interactions lncRNA FEELnc
#' @export
#' @importFrom utils read.table write.table
#' @examples
#' # --- 1. Create a temporary FEELnc classes file for a reproducible example ---
#' tempFile <- tempfile(fileext = ".tsv")
#'
#' # Create a sample data frame mimicking the FEELnc output structure
#' feelncData <- data.frame(
#'   lncRNA_transcript = c("lncT1", "lncT2", "lncT3", "lncT4"),
#'   lncRNA_gene = c("lncG1", "lncG1", "lncG2", "lncG2"),
#'   partnerRNA_transcript = c("mRNAT1", "mRNAT2", "mRNAT3", "mRNAT4"),
#'   partnerRNA_gene = c("mRNAG1", "mRNAG2", "mRNAG2", "mRNAG3"),
#'   distance = c(5000, 80000, 120000, 15000),
#'   isBest = c(1, 1, 0, 1)
#' )
#'
#' # Write the data to the temporary file
#' write.table(feelncData, file = tempFile, sep = "\t", quote = FALSE, row.names = FALSE)
#'
#' # --- 2. Run the cisInter function ---
#'
#' # Example 2a: Basic usage, filtering by 'isBest' and default distance
#' cisResults <- cisInter(feelncClasses = tempFile)
#' print("Results with default filtering:")
#' print(cisResults)
#'
#' # Example 2b: Advanced usage with ID and distance filtering
#' specificLncRNAs <- c("lncG1")
#' specificMRNAs <- c("mRNAG2")
#' cisResultsFiltered <- cisInter(
#'   feelncClasses = tempFile,
#'   lncRNAs = specificLncRNAs,
#'   lncRnaLevel = "gene",
#'   mRNAs = specificMRNAs,
#'   mRnaLevel = "gene",
#'   maxDist = 90000
#' )
#' print("Results with specific ID and distance filters:")
#' print(cisResultsFiltered)
#'
#' # --- 3. Clean up the temporary file ---
#' unlink(tempFile)
#'
cisInter <- function(feelncClasses, lncRNAs = NULL, mRNAs = NULL, isBest = TRUE, lncRnaLevel = "transcript", mRnaLevel = "gene", maxDist = 100000) {

  # --- 1. Input Validation ---
  if (missing(feelncClasses) || !is.character(feelncClasses) || length(feelncClasses) != 1) {
    stop("'feelncClasses' must be a single character string specifying a file path.")
  }
  if (!file.exists(feelncClasses)) {
    stop("The file specified in 'feelncClasses' does not exist: ", feelncClasses)
  }

  lncRnaLevel <- match.arg(lncRnaLevel, c("transcript", "gene"))
  mRnaLevel <- match.arg(mRnaLevel, c("transcript", "gene"))

  # --- 2. Read and Filter Data ---
  cis <- utils::read.table(feelncClasses, header = TRUE)

  if (isBest) {
    # Ensure the 'isBest' column exists before filtering
    if ("isBest" %in% colnames(cis)) {
      cis <- cis[cis$isBest == 1, ]
    } else {
      warning("'isBest' column not found in the input file. Filtering skipped.")
    }
  }

  if (!is.null(lncRNAs)) {
    if (lncRnaLevel == "transcript") {
      cis <- cis[cis$lncRNA_transcript %in% lncRNAs, ]
    } else {
      cis <- cis[cis$lncRNA_gene %in% lncRNAs, ]
    }
  }

  if (!is.null(mRNAs)) {
    if (mRnaLevel == "transcript") {
      cis <- cis[cis$partnerRNA_transcript %in% mRNAs, ]
    } else {
      cis <- cis[cis$partnerRNA_gene %in% mRNAs, ]
    }
  }

  cis <- cis[cis$distance <= maxDist, ]

  # --- 3. Rename Columns and Return ---
  # Only rename if the column exists to avoid errors
  if ("lncRNA_gene" %in% names(cis)) {
    names(cis)[names(cis) == 'lncRNA_gene'] <- 'lncRNAId'
  }
  if ("partnerRNA_gene" %in% names(cis)) {
    names(cis)[names(cis) == 'partnerRNA_gene'] <- 'targetRNAId'
  }

  return(cis)
}
