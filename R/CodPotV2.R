#' Combine Coding Potential Tool Results
#'
#' Reads outputs from various coding potential prediction tools and aggregates
#' the results into a list. Noncoding predictions are marked as 1, coding as 0.
#' Tool results are stored under a 'tools' sublist.
#'
#' @param CPC2_outfile Path to the CPC2 output file.
#' @param PLEK_outfile Path to the PLEK output file.
#' @param FEELnc_outfile Path to the FEELnc output file.
#' @param CPAT_outfile Path to the CPAT output file.
#' @param CPAT_cutoff Cutoff value for CPAT probability (default: 0.364).
#' @param CNCI_outfile Path to the CNCI output file.
#' @param LncFinder_outfile Path to the LncFinder output file.
#' @param lncRNA_Mdeep_outfile Path to the lncRNA_Mdeep output file.
#'
#' @return A list containing:
#'   \item{seqIDs}{A character vector of all unique sequence IDs found.}
#'   \item{tools}{A list where each element corresponds to a tool and contains
#'                a vector indicating noncoding (1) or coding (0) status for
#'                each sequence ID.}
#' @keywords coding potential lncRNA bioinformatics
#' @export
#' @examples
#' # Example usage (assuming output files exist):
#' # codpot_results <- CodPot2tbl(CPC2_outfile = "cpc2_results.txt",
#' #                              PLEK_outfile = "plek_results.txt")
#' # print(head(codpot_results$seqIDs))
#' # print(head(codpot_results$tools$CPC2))
#'
CodPot2tbl <- function(CPC2_outfile = NULL, PLEK_outfile = NULL, FEELnc_outfile = NULL, CPAT_outfile = NULL, CPAT_cutoff = 0.364, CNCI_outfile = NULL, LncFinder_outfile = NULL, lncRNA_Mdeep_outfile = NULL) {
  # Initialize a list to store results, with a dedicated sublist for tools
  codpot_list <- list(tools = list())

  # Read output files and store noncoding transcript IDs under the 'tools' sublist
  if (!is.null(CPC2_outfile)) {
    codpot_list$tools$CPC2 <- read.CPC2(CPC2_outfile)
  }
  if (!is.null(PLEK_outfile)) {
    codpot_list$tools$PLEK <- read.PLEK(PLEK_outfile)
  }
  if (!is.null(FEELnc_outfile)) {
    codpot_list$tools$FEELnc <- read.FEELnc(FEELnc_outfile)
  }
  if (!is.null(CPAT_outfile)) {
    codpot_list$tools$CPAT <- read.CPAT(CPAT_outfile, CPAT_cutoff)
  }
  if (!is.null(CNCI_outfile)) {
    codpot_list$tools$CNCI <- read.CNCI(CNCI_outfile)
  }
  if (!is.null(LncFinder_outfile)) {
    codpot_list$tools$LncFinder <- read.LncFinder(LncFinder_outfile)
  }
  if (!is.null(lncRNA_Mdeep_outfile)) {
    codpot_list$tools$lncRNA_Mdeep <- read.lncRNA_Mdeep(lncRNA_Mdeep_outfile)
  }

  # Get all unique sequence IDs from the tool results
  # unlist() combines all vectors within the 'tools' list into one
  seq_ids <- unique(unlist(codpot_list$tools))
  seq_ids <- seq_ids[seq_ids != ""]  # Remove empty strings
  seq_ids <- seq_ids[!is.na(seq_ids)] # Remove NA values

  # If no valid sequence IDs were found, return an empty list or handle appropriately
  if (length(seq_ids) == 0) {
    warning("No valid sequence IDs found in the provided files.")
    return(list(seqIDs = character(0), tools = list()))
  }

  # Create the final result structure: seqIDs and a matrix-like list for tool predictions
  result_list <- list(seqIDs = seq_ids)
  tool_predictions <- list()

  # Iterate through the tools and mark the presence (noncoding = 1) or absence (coding = 0) of each sequence ID
  for (tool_name in names(codpot_list$tools)) {
    tool_ids <- codpot_list$tools[[tool_name]]
    # ifelse checks if each seq_id is present in the current tool's noncoding list
    tool_predictions[[tool_name]] <- ifelse(seq_ids %in% tool_ids, 1, 0)
  }

  # Combine tool predictions into a matrix to easily check for all-zero columns
  if (length(tool_predictions) > 0) {
    prediction_matrix <- do.call(cbind, tool_predictions)
    # Keep only columns (tools) that predicted at least one noncoding transcript
    cols_to_keep <- colSums(prediction_matrix) > 0
    result_list$tools <- tool_predictions[cols_to_keep]
  } else {
    # Handle case where no tool files were provided or processed
    result_list$tools <- list()
  }


  return(result_list)
}


#' Read CPC2 Output
#'
#' Parses a CPC2 output file and extracts IDs of transcripts classified as noncoding.
#'
#' @param CPC2_outfile Path to the CPC2 output file.
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords CPC2 lncRNA parsing
#' @export
#' @examples
#' # noncoding_cpc2 <- read.CPC2("cpc2_results.txt")
#'
read.CPC2 <- function(CPC2_outfile) {
  # Read the table, assuming space-separated values and no header
  CPC2_data <- read.table(CPC2_outfile, stringsAsFactors = FALSE)
  # Filter rows where the 8th column indicates "noncoding"
  # Extract the transcript IDs from the 1st column for these rows
  CPC2_noncoding_ids <- CPC2_data[CPC2_data$V8 %in% "noncoding", ]$V1
  # Ensure IDs are unique and return as character vector
  CPC2_noncoding_ids <- unique(as.character(CPC2_noncoding_ids))
  return(CPC2_noncoding_ids)
}

#' Read PLEK Output
#'
#' Parses a PLEK output file and extracts IDs of transcripts classified as Non-coding.
#'
#' @param PLEK_outfile Path to the PLEK output file.
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords PLEK lncRNA parsing
#' @export
#' @examples
#' # noncoding_plek <- read.PLEK("plek_results.txt")
#'
read.PLEK <- function(PLEK_outfile) {
  # Read the table, assuming space-separated values and no header
  PLEK_data <- read.table(PLEK_outfile, stringsAsFactors = FALSE)
  # Remove the leading '>' character from the sequence IDs in the 3rd column
  PLEK_data$V3 <- gsub(">", "", PLEK_data$V3)
  # Filter rows where the 1st column indicates "Non-coding"
  # Extract the cleaned transcript IDs from the 3rd column for these rows
  PLEK_noncoding_ids <- PLEK_data[PLEK_data$V1 %in% "Non-coding", ]$V3
  # Ensure IDs are unique and return as character vector
  PLEK_noncoding_ids <- unique(as.character(PLEK_noncoding_ids))
  return(PLEK_noncoding_ids)
}

#' Read FEELnc Output
#'
#' Parses a FEELnc output file and extracts IDs of transcripts classified as noncoding (label 0).
#'
#' @param FEELnc_outfile Path to the FEELnc output file (table with 'label' column).
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords FEELnc lncRNA parsing
#' @export
#' @examples
#' # noncoding_feelnc <- read.FEELnc("feelnc_results.tab")
#'
read.FEELnc <- function(FEELnc_outfile) {
  # Read the table, assuming tab-separated values and presence of a header
  FEELnc_data <- read.table(FEELnc_outfile, header = TRUE, stringsAsFactors = FALSE)
  # Filter rows where the 'label' column is 0 (indicating noncoding)
  # Extract the transcript IDs from the 'name' column for these rows
  FEELnc_noncoding_ids <- FEELnc_data[FEELnc_data$label %in% 0, ]$name
  # Ensure IDs are unique and return as character vector
  FEELnc_noncoding_ids <- unique(as.character(FEELnc_noncoding_ids))
  return(FEELnc_noncoding_ids)
}

#' Read CPAT Output
#'
#' Parses a CPAT output file and extracts IDs of transcripts with a coding probability
#' below the specified cutoff.
#'
#' @param CPAT_outfile Path to the CPAT output file (table with 'coding_prob' column).
#' @param CPAT_cutoff The probability threshold below which transcripts are considered noncoding.
#'                    Default is 0.364 (predefined for human).
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords CPAT lncRNA parsing cutoff
#' @export
#' @examples
#' # noncoding_cpat <- read.CPAT("cpat_results.txt", CPAT_cutoff = 0.4)
#'
read.CPAT <- function(CPAT_outfile, CPAT_cutoff = 0.364) {
  # Read the table, assuming space/tab-separated values and presence of a header
  # row.names=1 might be needed if IDs are row names, adjust if IDs are in a column
  CPAT_data <- read.table(CPAT_outfile, header = TRUE, stringsAsFactors = FALSE) # Potentially add row.names=1 if needed
  # Filter rows where 'coding_prob' is less than the cutoff
  # Extract the transcript IDs (assuming they are row names)
  # If IDs are in a specific column (e.g., "transcript_id"), use CPAT_data[..., "transcript_id"]
  CPAT_noncoding_ids <- rownames(CPAT_data[CPAT_data$coding_prob < CPAT_cutoff, ])
  # Ensure IDs are unique and return as character vector (already character from rownames)
  CPAT_noncoding_ids <- unique(CPAT_noncoding_ids)
  return(CPAT_noncoding_ids)
}

#' Read CNCI Output
#'
#' Parses a CNCI output file and extracts IDs of transcripts classified as noncoding.
#'
#' @param CNCI_outfile Path to the CNCI output file.
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords CNCI lncRNA parsing
#' @export
#' @examples
#' # noncoding_cnci <- read.CNCI("cnci_results.tsv")
#'
read.CNCI <- function(CNCI_outfile) {
  # Read the table, assuming tab-separated values and presence of a header
  CNCI_data <- read.table(CNCI_outfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Filter rows where the 'index' column indicates "noncoding"
  # Extract the transcript IDs from the first column for these rows
  CNCI_noncoding_ids <- CNCI_data[CNCI_data$index %in% "noncoding", 1]
  # Ensure IDs are unique and return as character vector
  CNCI_noncoding_ids <- unique(as.character(CNCI_noncoding_ids))
  return(CNCI_noncoding_ids)
}

#' Read LncFinder Output
#'
#' Parses a LncFinder output file (CSV format) and extracts IDs of transcripts
#' predicted as "NonCoding".
#'
#' @param LncFinder_outfile Path to the LncFinder output file (CSV format).
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords LncFinder lncRNA parsing csv
#' @export
#' @examples
#' # noncoding_lncfinder <- read.LncFinder("lncfinder_results.csv")
#'
read.LncFinder <- function(LncFinder_outfile) {
  # Read the CSV file, assuming semicolon separator, header, and IDs as row names (adjust if needed)
  LncFinder_data <- read.csv2(LncFinder_outfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  # Filter rows where the 'Pred' column indicates "NonCoding"
  # Extract the transcript IDs (assuming they are row names)
  LncFinder_noncoding_ids <- rownames(LncFinder_data[LncFinder_data$Pred %in% "NonCoding", ])
  # Ensure IDs are unique and return as character vector (already character from rownames)
  LncFinder_noncoding_ids <- unique(LncFinder_noncoding_ids)
  return(LncFinder_noncoding_ids)
}

#' Read lncRNA_Mdeep Output
#'
#' Parses an lncRNA_Mdeep output file and extracts IDs of transcripts classified as noncoding.
#'
#' @param lncRNA_Mdeep_outfile Path to the lncRNA_Mdeep output file.
#'
#' @return A character vector of noncoding transcript IDs.
#' @keywords lncRNA_Mdeep lncRNA parsing
#' @export
#' @examples
#' # noncoding_mdeep <- read.lncRNA_Mdeep("mdeep_results.txt")
#'
read.lncRNA_Mdeep <- function(lncRNA_Mdeep_outfile) {
  # Read the table, assuming space-separated values and no header
  lncRNA_Mdeep_data <- read.table(lncRNA_Mdeep_outfile, stringsAsFactors = FALSE)
  # Filter rows where the second column indicates "noncoding"
  # Extract the transcript IDs from the first column for these rows
  lncRNA_Mdeep_noncoding_ids <- lncRNA_Mdeep_data[lncRNA_Mdeep_data$V2 %in% "noncoding", ]$V1
  # Ensure IDs are unique and return as character vector
  lncRNA_Mdeep_noncoding_ids <- unique(as.character(lncRNA_Mdeep_noncoding_ids))
  return(lncRNA_Mdeep_noncoding_ids)
}


#' Create Venn Diagram from Coding Potential Results
#'
#' Generates a Venn diagram illustrating the overlap of noncoding predictions
#' from different tools, based on the output of `CodPot2tbl`.
#'
#' @param CodPot A list object generated by `CodPot2tbl()`. Must contain
#'               `$seqIDs` and `$tools` sublist.
#' @param selmet An optional numeric vector indicating which tools to include
#'               in the diagram (1 for include, 0 for exclude). Order corresponds
#'               to the order of tools in `CodPot$tools`. If NULL, all tools are included.
#' @param venncolors A vector of colors for the Venn diagram segments.
#'
#' @return Invisibly returns the result of `venn::venn`, which might be used
#'         for further plot manipulation if needed. The primary output is the plot.
#' @keywords venn diagram coding potential visualization plot lncRNA
#' @export
#' @importFrom venn venn
#' @examples
#' # Assuming 'codpot_results' is output from CodPot2tbl() with CPC2, PLEK, CPAT
#' # venn.CodPot(CodPot = codpot_results) # Venn for all tools
#' # venn.CodPot(CodPot = codpot_results, selmet = c(1, 0, 1)) # Venn for CPC2 and CPAT only
#'
venn.CodPot <- function(CodPot, venncolors = c("green", "red", "blue", "yellow", "magenta", "black", "orange", "white", "darkgreen"), selmet = NULL) {
  # Check if the input structure is as expected
  if (!is.list(CodPot) || !("seqIDs" %in% names(CodPot)) || !("tools" %in% names(CodPot)) || !is.list(CodPot$tools)) {
    stop("Input 'CodPot' must be a list generated by CodPot2tbl() with 'seqIDs' and 'tools' elements.")
  }

  num_tools <- length(CodPot$tools)
  if (num_tools == 0) {
    warning("No tool results found in 'CodPot$tools'. Cannot generate Venn diagram.")
    return(invisible(NULL))
  }

  # If selmet is not provided, include all tools
  if (is.null(selmet)) {
    selmet <- rep(1, num_tools)
  } else if (length(selmet) != num_tools) {
    stop(paste("Length of 'selmet' (", length(selmet), ") does not match the number of tools (", num_tools, ")."))
  }

  # Prepare the list for the venn function based on selected methods
  lista_for_venn <- list()
  selected_tool_names <- names(CodPot$tools)[which(selmet == 1)]

  # Check if any tools were selected
  if (length(selected_tool_names) == 0) {
    warning("No tools selected by 'selmet'. Cannot generate Venn diagram.")
    return(invisible(NULL))
  }

  # Populate the list with the actual sequence IDs for selected tools
  for (tool_name in selected_tool_names) {
    # Get the binary vector (0/1) for the tool
    tool_binary_vector <- CodPot$tools[[tool_name]]
    # Extract the seqIDs where the value is 1 (noncoding)
    lista_for_venn[[tool_name]] <- CodPot$seqIDs[which(tool_binary_vector == 1)]
  }

  # Generate the Venn diagram using the 'venn' package
  # Ensure enough colors are available for the selected tools
  num_selected <- sum(selmet)
  if (length(venncolors) < num_selected) {
    warning("Not enough colors provided in 'venncolors'. Recycling colors.")
    # Simple recycling, could be improved with better color generation if needed
    venncolors <- rep(venncolors, length.out = num_selected)
  }

  # Call the venn function
  venn_result <- venn::venn(
    lista_for_venn,
    zcolor = venncolors[1:num_selected],
    snames = selected_tool_names, # Use the names of the selected tools
    ilabels = TRUE # Show intersection counts
  )

  # Return the result invisibly (the plot is the main output)
  invisible(venn_result)
}
