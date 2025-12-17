#' Find Potential cis-Regulatory Interactions from FEELnc Output
#'
#' This function reads and filters the output file from FEELnc (`.classes` file)
#' to identify potential cis-interactions between lncRNAs and mRNAs based on
#' user-specified criteria such as genomic distance and specific gene/transcript lists.
#'
#' @param FEELncClassesFile A character string specifying the path to the
#'   FEELnc `.classes` output file. This file is required.
#' @param lncRnas An optional character vector of lncRNA IDs (gene or transcript)
#'   to filter the results. If `NULL`, no filtering by lncRNA ID is performed.
#' @param mRnas An optional character vector of mRNA IDs (gene or transcript)
#'   to filter the results. If `NULL`, no filtering by mRNA ID is performed.
#' @param filterIsBest Logical. If `TRUE` (default), only interactions marked as
#'   "best" (`isBest == 1`) in the FEELnc output are kept.
#' @param lncRnaLevel A character string specifying the level for lncRNA filtering:
#'   `"transcript"` (default) or `"gene"`.
#' @param mRnaLevel A character string specifying the level for mRNA filtering:
#'   `"gene"` (default) or `"transcript"`.
#' @param maxDist A numeric value for the maximum distance (in base pairs) to
#'   consider for a cis-interaction (default: 100,000).
#'
#' @return A `data.frame` containing the filtered cis-interaction data with
#'   renamed columns (`lncRNAId`, `targetRNAId`).
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' # --- 1. Create a temporary FEELnc output file for a reproducible example ---
#' feelncFile <- tempfile()
#' mock_data <- data.frame(
#'   isBest = c(1, 1, 0, 1, 1),
#'   lncRNA_gene = c("LNC_G1", "LNC_G1", "LNC_G2", "LNC_G3", "LNC_G4"),
#'   lncRNA_transcript = c("LNC_T1", "LNC_T2", "LNC_T3", "LNC_T4", "LNC_T5"),
#'   partnerRNA_gene = c("TARGET_G1", "TARGET_G2", "TARGET_G1", "TARGET_G3", "TARGET_G4"),
#'   partnerRNA_transcript = c("T_T1", "T_T2", "T_T3", "T_T4", "T_T5"),
#'   distance = c(5000, 80000, 1000, 120000, 9000)
#' )
#' utils::write.table(mock_data, feelncFile, sep = "\t", row.names = FALSE, col.names = TRUE)
#'
#' # --- 2. Define lncRNA and mRNA lists for filtering ---
#' lncRnaList <- c("LNC_T1", "LNC_T4") # Filter by transcript ID
#' mRnaList <- c("TARGET_G1")       # Filter by gene ID
#'
#' # --- 3. Run the function ---
#' cis_interactions <- findCisInteractions(
#'   FEELncClassesFile = feelncFile,
#'   lncRnas = lncRnaList,
#'   mRnas = mRnaList,
#'   lncRnaLevel = "transcript",
#'   mRnaLevel = "gene",
#'   maxDist = 100000
#' )
#'
#' print(cis_interactions)
#'
#' # --- 4. Clean up the temporary file ---
#' unlink(feelncFile)
#'
findCisInteractions <- function(FEELncClassesFile, lncRnas = NULL, mRnas = NULL,
                                filterIsBest = TRUE, lncRnaLevel = "transcript",
                                mRnaLevel = "gene", maxDist = 100000) {

    # --- 1. Input Validation ---
    if (!file.exists(FEELncClassesFile)) {
        stop("The specified 'FEELncClassesFile' does not exist: ", FEELncClassesFile)
    }
    lncRnaLevel <- match.arg(lncRnaLevel, choices = c("transcript", "gene"))
    mRnaLevel <- match.arg(mRnaLevel, choices = c("transcript", "gene"))

    # --- 2. Read and Pre-filter Data ---
    cis <- utils::read.table(FEELncClassesFile, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)

    if (filterIsBest) {
        if ("isBest" %in% colnames(cis)) {
            cis <- cis[cis$isBest == 1, ]
        } else {
            warning("'isBest' column not found, skipping this filter.", call. = FALSE)
        }
    }

    # --- 3. Filter by lncRNA and mRNA lists ---
    if (!is.null(lncRnas) && length(lncRnas) > 0) {
        col_to_filter <- if (lncRnaLevel == "transcript") "lncRNA_transcript" else "lncRNA_gene"
        if (col_to_filter %in% colnames(cis)) {
            cis <- cis[cis[[col_to_filter]] %in% lncRnas, ]
        }
    }

    if (!is.null(mRnas) && length(mRnas) > 0) {
        col_to_filter <- if (mRnaLevel == "transcript") "partnerRNA_transcript" else "partnerRNA_gene"
        if (col_to_filter %in% colnames(cis)) {
            cis <- cis[cis[[col_to_filter]] %in% mRnas, ]
        }
    }

    # --- 4. Filter by Distance and Rename Columns ---
    if ("distance" %in% colnames(cis)) {
        cis <- cis[cis$distance <= maxDist, ]
    }

    # Rename columns to a consistent style
    current_names <- colnames(cis)
    current_names[current_names == 'lncRNA_gene'] <- 'lncRNAId'
    current_names[current_names == 'partnerRNA_gene'] <- 'targetRNAId'
    colnames(cis) <- current_names

    return(cis)
}
