#' Run LION RNA-Protein Interaction Analysis
#'
#' Performs RNA-protein interaction predictions using the LION package.
#' This function takes lists of RNA and protein sequences and runs the analysis
#' on all possible pairs, saving the results to a CSV file.
#'
#' @param rnaSeqs A named list of RNA sequences, where names are sequence IDs
#'   and values are the sequences themselves (as strings).
#' @param protSeqs A named list of protein sequences, structured similarly to
#'   `rnaSeqs`.
#' @param outputFilename A character string specifying the path for the output
#'   CSV file. Defaults to "LION_results.csv".
#' @param parallelCores An integer specifying the number of CPU cores for
#'   parallel processing. Defaults to 2.
#'
#' @return A data.frame containing the LION prediction results for all
#'   RNA-protein pairs.
#'
#' @importFrom LION run_confidentPrediction
#' @importFrom utils write.csv2
#' @export
#' @examples
#' # This is a time-consuming function. The example demonstrates its usage
#' # by "mocking" the core LION function to run quickly during checks.
#'
#' # --- 1. Create sample sequence data ---
#' sampleRnaSeqs <- list(
#'   RNA_1 = "AUGGCUAGU",
#'   RNA_2 = "GCUAGUAGC"
#' )
#' sampleProtSeqs <- list(
#'   PROT_A = "MKTAY",
#'   PROT_B = "LFWDP"
#' )
#'
#' # --- 2. Mock the time-consuming function for this example ---
#' # In a real analysis, the actual LION::run_confidentPrediction would be called.
#' if (requireNamespace("LION", quietly = TRUE)) {
#'   # Define a fake function that returns a plausible result instantly.
#'   run_confidentPrediction <- function(seqRNA, seqPro, methods, label, parallel.cores) {
#'     # Generate all pairs from the input lists
#'     pairs <- expand.grid(RNA = names(seqRNA), Protein = names(seqPro))
#'     # Create a dummy result data.frame
#'     mock_result <- data.frame(
#'       RNA = pairs$RNA,
#'       Protein = pairs$Protein,
#'       LION_Score = round(runif(nrow(pairs)), 3)
#'     )
#'     # LION returns a list, so we mimic that structure
#'     return(list(mock_result))
#'   }
#'
#'   # --- 3. Run the analysis with the mocked function ---
#'   tempFile <- tempfile(fileext = ".csv")
#'   results <- runLionAnalysis(
#'     rnaSeqs = sampleRnaSeqs,
#'     protSeqs = sampleProtSeqs,
#'     outputFilename = tempFile,
#'     parallelCores = 1 # Use 1 core for the example
#'   )
#'
#'   print("Mock analysis results:")
#'   print(results)
#'
#'   # --- 4. Clean up ---
#'   unlink(tempFile)
#' }
runLionAnalysis <- function(rnaSeqs, protSeqs,
                            outputFilename = "LION_results.csv",
                            parallelCores = 2) {

  # --- 1. Input Validation ---
  stopifnot(
    "'rnaSeqs' must be a named list with at least one sequence" =
      (is.list(rnaSeqs) && !is.null(names(rnaSeqs)) && length(rnaSeqs) > 0),
    "'protSeqs' must be a named list with at least one sequence" =
      (is.list(protSeqs) && !is.null(names(protSeqs)) && length(protSeqs) > 0)
  )

  # --- 2. Run LION Prediction (Efficiently) ---
  # The core LION function is designed to handle all pairs internally.
  # The original looping logic was inefficient and has been removed.
  message("Starting LION RNA-Protein Interaction Analysis...")
  message("Analyzing ", length(rnaSeqs), " RNAs against ", length(protSeqs), " proteins.")

  predictionResults <- LION::run_confidentPrediction(
    seqRNA = rnaSeqs,
    seqPro = protSeqs,
    methods = c("RPISeq_retrain", "rpiCOOL_retrain", "LION"),
    label = "Interact",
    parallel.cores = parallelCores
  )

  # --- 3. Process and Save Results ---
  # Combine results into a single data.frame and remove duplicated columns
  predictionDf <- as.data.frame(do.call("cbind", predictionResults))
  predictionDf <- predictionDf[, !duplicated(colnames(predictionDf))]

  utils::write.csv2(predictionDf, outputFilename, row.names = FALSE)
  message("Analysis completed. Results saved to: ", outputFilename)

  return(predictionDf)
}
