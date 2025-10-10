#' Combine Coding Potential Tool Results
#'
#' Reads outputs from various tools and aggregates results into a structured
#' list. Noncoding predictions are marked as 1, coding as 0.
#'
#' @param cpc2Outfile Path to the CPC2 output file.
#' @param plekOutfile Path to the PLEK output file.
#' @param feelncOutfile Path to the FEELnc output file.
#' @param cpatOutfile Path to the CPAT output file.
#' @param cpatCutoff Cutoff value for CPAT probability (default: 0.364).
#' @param cnciOutfile Path to the CNCI output file.
#' @param lncFinderOutfile Path to the LncFinder output file.
#' @param lncRnaMdeepOutfile Path to the lncRNA_Mdeep output file.
#'
#' @return A list containing `seqIDs` (a character vector of all unique
#'   sequence IDs) and `tools` (a list of binary vectors for each tool).
#' @keywords coding potential lncRNA bioinformatics
#' @export
#' @importFrom utils read.table write.table
#' @examples
#' # --- 1. Create temporary output files for a reproducible example ---
#' cpc2File <- tempfile()
#' plekFile <- tempfile()
#' # CPC2 format: ID ... classification
#' write.table(
#'   data.frame(V1="ID1", V2="-", V3="-", V4="-", V5="-", V6="-", V7="-", V8="noncoding"),
#'   cpc2File, col.names=FALSE, row.names=FALSE, sep="\t"
#' )
#' write.table(
#'   data.frame(V1="ID2", V2="-", V3="-", V4="-", V5="-", V6="-", V7="-", V8="coding"),
#'   cpc2File, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE
#' )
#' # PLEK format: Classification Score >ID
#' write.table(
#'   data.frame(V1="Non-coding", V2="0.9", V3=">ID2"),
#'   plekFile, col.names=FALSE, row.names=FALSE, sep=" "
#' )
#'
#' # --- 2. Run the function with the temporary files ---
#' codpotResults <- codPotToTbl(cpc2Outfile = cpc2File, plekOutfile = plekFile)
#'
#' print("Sequence IDs:")
#' print(codpotResults$seqIDs)
#' print("Tool predictions:")
#' print(codpotResults$tools)
#'
#' # --- 3. Clean up the temporary files ---
#' unlink(c(cpc2File, plekFile))
#'
codPotToTbl <- function(cpc2Outfile = NULL, plekOutfile = NULL, feelncOutfile = NULL,
                        cpatOutfile = NULL, cpatCutoff = 0.364, cnciOutfile = NULL,
                        lncFinderOutfile = NULL, lncRnaMdeepOutfile = NULL) {

  codpotList <- list(tools = list())

  # Helper to ensure file exists before reading
  readFileIfExists <- function(outfile, readFunc, ...) {
    if (!is.null(outfile)) {
      if (file.exists(outfile)) {
        return(readFunc(outfile, ...))
      } else {
        warning("File not found: ", outfile, ". Skipping.", call. = FALSE)
        return(NULL)
      }
    }
    return(NULL)
  }

  codpotList$tools$CPC2 <- readFileIfExists(cpc2Outfile, readCpc2)
  codpotList$tools$PLEK <- readFileIfExists(plekOutfile, readPlek)
  codpotList$tools$FEELnc <- readFileIfExists(feelncOutfile, readFeelnc)
  codpotList$tools$CPAT <- readFileIfExists(cpatOutfile, readCpat, CPAT_cutoff = cpatCutoff)
  codpotList$tools$CNCI <- readFileIfExists(cnciOutfile, readCnci)
  codpotList$tools$LncFinder <- readFileIfExists(lncFinderOutfile, readLncFinder)
  codpotList$tools$lncRNA_Mdeep <- readFileIfExists(lncRnaMdeepOutfile, readLncRnaMdeep)

  # Remove NULL entries from tools that were not provided or found
  codpotList$tools <- Filter(Negate(is.null), codpotList$tools)

  if (length(codpotList$tools) == 0) {
    warning("No valid tool output files were provided or found.")
    return(list(seqIDs = character(0), tools = list()))
  }

  seqIds <- unique(unlist(lapply(codpotList$tools, names)))
  seqIds <- seqIds[!is.na(seqIds) & seqIds != ""]

  if (length(seqIds) == 0) {
    warning("No valid sequence IDs found in the provided files.")
    return(list(seqIDs = character(0), tools = list()))
  }

  resultList <- list(seqIDs = seqIds)
  toolPredictions <- lapply(codpotList$tools, function(toolIds) {
    # The output from reader functions is now a named vector
    ifelse(seqIds %in% names(toolIds), 1, 0)
  })

  predictionMatrix <- do.call(cbind, toolPredictions)
  colsToKeep <- colSums(predictionMatrix, na.rm = TRUE) > 0
  resultList$tools <- toolPredictions[colsToKeep]

  return(resultList)
}


#' Read CPC2 Output
#' @param cpc2Outfile Path to the CPC2 output file.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @importFrom utils read.table write.table
#' @examples
#' # Create a temporary file for the example
#' tempCpc2File <- tempfile()
#' write.table(
#'   data.frame(V1=c("ID1", "ID2", "ID3"), V8=c("noncoding", "coding", "noncoding")),
#'   tempCpc2File, row.names=FALSE, col.names=FALSE, quote=FALSE
#' )
#' # Run the reader function
#' noncodingCpc2 <- readCpc2(tempCpc2File)
#' print(noncodingCpc2)
#' # Clean up
#' unlink(tempCpc2File)
readCpc2 <- function(cpc2Outfile) {
  cpc2Data <- utils::read.table(cpc2Outfile, stringsAsFactors = FALSE)
  noncodingIds <- cpc2Data[cpc2Data$V8 == "noncoding", 1]
  noncodingIds <- unique(as.character(noncodingIds))
  # Return a named vector for easier processing later
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read PLEK Output
#' @param plekOutfile Path to the PLEK output file.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary file
#' tempPlekFile <- tempfile()
#' write.table(
#'  data.frame(V1=c("Non-coding", "Coding"), V3=c(">ID_A", ">ID_B")),
#'  tempPlekFile, row.names=FALSE, col.names=FALSE, quote=FALSE
#' )
#' # Run and print
#' noncodingPlek <- readPlek(tempPlekFile)
#' print(noncodingPlek)
#' # Clean up
#' unlink(tempPlekFile)
readPlek <- function(plekOutfile) {
  plekData <- utils::read.table(plekOutfile, header = FALSE, stringsAsFactors = FALSE)

  if (ncol(plekData) < 3) {
    warning("PLEK output file '", plekOutfile, "' has fewer than 3 columns. Cannot extract IDs.")
    return(setNames(character(0), character(0)))
  }

  plekData$V3 <- gsub(">", "", plekData$V3)
  noncodingIds <- plekData[plekData$V1 == "Non-coding", 3]
  noncodingIds <- unique(as.character(noncodingIds))
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read FEELnc Output
#' @param feelncOutfile Path to the FEELnc output file.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary file
#' tempFeelncFile <- tempfile()
#' write.table(
#'  data.frame(name=c("T1", "T2", "T3"), label=c(0, 1, 0)),
#'  tempFeelncFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t"
#' )
#' # Run and print
#' noncodingFeelnc <- readFeelnc(tempFeelncFile)
#' print(noncodingFeelnc)
#' # Clean up
#' unlink(tempFeelncFile)
readFeelnc <- function(feelncOutfile) {
  feelncData <- utils::read.table(feelncOutfile, header = TRUE, stringsAsFactors = FALSE)
  noncodingIds <- feelncData[feelncData$label == 0, "name"]
  noncodingIds <- unique(as.character(noncodingIds))
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read CPAT Output
#' @param cpatOutfile Path to the CPAT output file.
#' @param cpatCutoff Probability threshold for noncoding classification.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary file
#' tempCpatFile <- tempfile()
#' cpat_df <- data.frame(coding_prob=c(0.1, 0.9, 0.3))
#' rownames(cpat_df) <- c("TX1", "TX2", "TX3")
#' write.table(cpat_df, tempCpatFile, quote=FALSE)
#' # Run and print
#' noncodingCpat <- readCpat(tempCpatFile, cpatCutoff = 0.4)
#' print(noncodingCpat)
#' # Clean up
#' unlink(tempCpatFile)
readCpat <- function(cpatOutfile, cpatCutoff = 0.364) {
  cpatData <- utils::read.table(cpatOutfile, header = TRUE, stringsAsFactors = FALSE)
  noncodingIds <- rownames(cpatData[cpatData$coding_prob < cpatCutoff, ])
  noncodingIds <- unique(noncodingIds)
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read CNCI Output
#' @param cnciOutfile Path to the CNCI output file.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary file
#' tempCnciFile <- tempfile()
#' write.table(
#'   data.frame(transcript_ID=c("id_x", "id_y"), index=c("noncoding", "coding")),
#'   tempCnciFile, row.names=FALSE, quote=FALSE, sep="\t"
#' )
#' # Run and print
#' noncodingCnci <- readCnci(tempCnciFile)
#' print(noncodingCnci)
#' # Clean up
#' unlink(tempCnciFile)
readCnci <- function(cnciOutfile) {
  cnciData <- utils::read.table(cnciOutfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  noncodingIds <- cnciData[cnciData$index == "noncoding", 1]
  noncodingIds <- unique(as.character(noncodingIds))
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read LncFinder Output
#' @param lncFinderOutfile Path to the LncFinder output file (CSV).
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @importFrom utils read.csv2
#' @examples
#' # Create a temporary file
#' tempLncFinderFile <- tempfile()
#' lncf_df <- data.frame(ID=c("rna1", "rna2"), Pred=c("NonCoding", "Coding"))
#' # write.csv2 is semicolon separated
#' write.csv2(lncf_df, tempLncFinderFile, row.names=FALSE, quote=FALSE)
#' # Adjust reader to expect 'ID' as first column instead of rownames
#' # Run and print
#' noncodingLncFinder <- readLncFinder(tempLncFinderFile)
#' print(noncodingLncFinder)
#' # Clean up
#' unlink(tempLncFinderFile)
readLncFinder <- function(lncFinderOutfile) {
  # LncFinder often uses the first column for IDs, not rownames
  lncFinderData <- utils::read.csv2(lncFinderOutfile, header = TRUE, stringsAsFactors = FALSE)
  # Assuming the first column contains the IDs
  noncodingIds <- lncFinderData[lncFinderData$Pred == "NonCoding", 1]
  noncodingIds <- unique(as.character(noncodingIds))
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' Read lncRNA_Mdeep Output
#' @param lncRnaMdeepOutfile Path to the lncRNA_Mdeep output file.
#' @return A named character vector of noncoding transcript IDs.
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary file
#' tempMdeepFile <- tempfile()
#' write.table(
#'  data.frame(V1=c("seqA", "seqB"), V2=c("noncoding", "coding")),
#'  tempMdeepFile, row.names=FALSE, col.names=FALSE, quote=FALSE
#' )
#' # Run and print
#' noncodingMdeep <- readLncRnaMdeep(tempMdeepFile)
#' print(noncodingMdeep)
#' # Clean up
#' unlink(tempMdeepFile)
readLncRnaMdeep <- function(lncRnaMdeepOutfile) {
  lncRnaMdeepData <- utils::read.table(lncRnaMdeepOutfile, stringsAsFactors = FALSE)
  noncodingIds <- lncRnaMdeepData[lncRnaMdeepData$V2 == "noncoding", 1]
  noncodingIds <- unique(as.character(noncodingIds))
  setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}


#' Create Venn Diagram from Coding Potential Results
#'
#' Generates a Venn diagram of noncoding predictions from tools, based on the
#' output of `codPotToTbl`. This function requires the 'venn' package.
#'
#' @param codPot A list object generated by `codPotToTbl()`.
#' @param selection An optional numeric vector indicating which tools to include
#'   (1 for include, 0 for exclude). If `NULL`, all tools are included.
#' @param vennColors A vector of colors for the Venn diagram segments.
#'
#' @return Invisibly returns the result of `venn::venn`. The primary output is
#'   the generated plot. If the 'venn' package is not installed, it returns
#'   `invisible(NULL)` after printing a message.
#' @keywords venn diagram coding potential visualization plot lncRNA
#' @export
#' @examples
#' # --- 1. Create a mock object that mimics the output of codPotToTbl ---
#' mockCodPot <- list(
#'   seqIDs = c("ID1", "ID2", "ID3", "ID4"),
#'   tools = list(
#'     CPC2 = c(1, 1, 0, 0), # ID1, ID2
#'     PLEK = c(0, 1, 1, 0), # ID2, ID3
#'     CPAT = c(0, 0, 1, 1)  # ID3, ID4
#'   )
#' )
#'
#' # --- 2. Generate Venn diagrams ---
#' # The example will only run if the 'venn' package is installed.
#' if (requireNamespace("venn", quietly = TRUE)) {
#'   # To prevent plots from showing up during automated checks, we can
#'   # temporarily redirect the graphics output.
#'   png(tempfile())
#'
#'   # Example with all tools
#'   plotCodPotVenn(codPot = mockCodPot)
#'
#'   # Example with a selection of tools (CPC2 and CPAT)
#'   plotCodPotVenn(codPot = mockCodPot, selection = c(1, 0, 1))
#'
#'   dev.off() # Close the graphics device
#' } else {
#'   message("Skipping venn diagram examples because the 'venn' package is not installed.")
#' }
plotCodPotVenn <- function(codPot, selection = NULL,
                           vennColors = c("green", "red", "blue", "yellow", "magenta",
                                          "black", "orange", "white", "darkgreen")) {

  # --- Check for suggested package ---
  if (!requireNamespace("venn", quietly = TRUE)) {
    message("The 'venn' package is required for the plotCodPotVenn() function.")
    message("Please install it using: BiocManager::install('venn')")
    return(invisible(NULL)) # Gracefully exit if the package is missing
  }

  stopifnot(
    "Input 'codPot' must be a list" = is.list(codPot),
    "'codPot' must contain 'seqIDs' and 'tools'" = all(c("seqIDs", "tools") %in% names(codPot)),
    "'codPot$tools' must be a list" = is.list(codPot$tools)
  )

  numTools <- length(codPot$tools)
  if (numTools == 0) {
    warning("No tool results found. Cannot generate Venn diagram.")
    return(invisible(NULL))
  }

  if (is.null(selection)) {
    selection <- rep(1, numTools)
  } else if (length(selection) != numTools) {
    stop(paste("Length of 'selection' (", length(selection),
               ") does not match the number of tools (", numTools, ")."))
  }

  selectedToolNames <- names(codPot$tools)[which(selection == 1)]

  if (length(selectedToolNames) == 0) {
    warning("No tools selected. Cannot generate Venn diagram.")
    return(invisible(NULL))
  }

  listForVenn <- lapply(selectedToolNames, function(toolName) {
    toolBinaryVector <- codPot$tools[[toolName]]
    codPot$seqIDs[which(toolBinaryVector == 1)]
  })
  names(listForVenn) <- selectedToolNames

  numSelected <- sum(selection)
  if (length(vennColors) < numSelected) {
    warning("Not enough colors provided; recycling colors.")
    vennColors <- rep(vennColors, length.out = numSelected)
  }

  vennResult <- venn::venn(
    listForVenn,
    zcolor = vennColors[seq_len(numSelected)],
    snames = selectedToolNames,
    ilabels = TRUE
  )

  invisible(vennResult)
}
