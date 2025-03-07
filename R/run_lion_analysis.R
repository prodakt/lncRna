#' Run LION RNA-Protein Interaction Analysis
#'
#' Performs RNA-protein interaction predictions using the LION package by processing all
#' combinations of RNA and protein sequences in pairs. Results are saved iteratively to a CSV file.
#'
#' @param rna_seqs Named list of RNA sequences, typically from \code{seqinr::read.fasta}.
#' @param prot_seqs Named list of protein sequences, typically from \code{seqinr::read.fasta}.
#' @param output_filename Character string specifying the output CSV file name (default: "LION_results.csv").
#' @param parallel_cores Integer specifying the number of CPU cores for parallel processing (default: 8).
#' @return A data frame containing the cumulative LION prediction results.
#' @export
#' @examples
#' \dontrun{
#'   # Load example sequence files
#'   rna_seqs <- seqinr::read.fasta("rna_sequences.fa", seqtype = "RNA", as.string = TRUE)
#'   prot_seqs <- seqinr::read.fasta("protein_sequences.fa", seqtype = "AA", as.string = TRUE)
#'   
#'   # Run LION analysis
#'   results <- run_lion_analysis(rna_seqs, prot_seqs, output_filename = "my_results.csv", parallel_cores = 4)
#' }
run_lion_analysis <- function(rna_seqs, prot_seqs, output_filename = "LION_results.csv", parallel_cores = 8) {
  # Validate input sequence lengths
  if (length(rna_seqs) < 2 | length(prot_seqs) < 2) {
    stop("Each sequence input must contain at least 2 sequences for LION analysis!")
  }
  
  # Generate all RNA-protein combinations
  combinations <- expand.grid(RNA = names(rna_seqs), Protein = names(prot_seqs))
  
  # Ensure even number of combinations
  if (nrow(combinations) %% 2 == 1) {
    combinations <- rbind(combinations, combinations[nrow(combinations), ])
  }
  
  total_sets <- nrow(combinations) / 2
  set_counter <- 0
  lion_results <- NULL
  
  cat("Starting LION RNA-Protein Interaction Analysis...\n")
  cat("Total sets to analyze:", total_sets, "\n")
  
  # Process combinations in pairs
  for (i in seq(1, nrow(combinations), by = 2)) {
    rna_name_1 <- combinations$RNA[i]
    rna_name_2 <- combinations$RNA[i + 1]
    prot_name_1 <- combinations$Protein[i]
    prot_name_2 <- combinations$Protein[i + 1]
    
    rna_subset <- list(rna_seqs[[rna_name_1]], rna_seqs[[rna_name_2]])
    names(rna_subset) <- c(rna_name_1, rna_name_2)
    
    prot_subset <- list(prot_seqs[[prot_name_1]], prot_seqs[[prot_name_2]])
    names(prot_subset) <- c(prot_name_1, prot_name_2)
    
    set_counter <- set_counter + 1
    cat("\nAnalyzing set:", set_counter, "/", total_sets, "\n")
    cat(sprintf("Progress: %.2f%%\n", (set_counter / total_sets) * 100))
    cat("RNA:", rna_name_1, "vs. Protein:", prot_name_1, "\n")
    cat("RNA:", rna_name_2, "vs. Protein:", prot_name_2, "\n")
    
    # Run LION prediction
    prediction_results <- run_confidentPrediction(
      seqRNA = rna_subset,
      seqPro = prot_subset,
      methods = c("RPISeq_retrain", "rpiCOOL_retrain", "LION"),
      label = "Interact",
      parallel.cores = parallel_cores
    )
    
    # Convert to data frame and deduplicate columns
    prediction_df <- do.call("cbind", prediction_results)
    prediction_df <- prediction_df[, !duplicated(colnames(prediction_df))]
    
    # Accumulate results
    lion_results <- if (is.null(lion_results)) prediction_df else rbind(lion_results, prediction_df)
    
    # Save results incrementally
    write.csv2(lion_results, output_filename, row.names = FALSE)
    cat("Results saved to", output_filename, "\n")
  }
  
  cat("\nAnalysis completed. Results saved to:", getwd(), "/", output_filename, "\n")
  return(lion_results)
}