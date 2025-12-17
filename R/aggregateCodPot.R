#' Aggregate Coding Potential Tool Results
#'
#' Reads output files from various coding potential prediction tools and
#' aggregates their results into a structured list. Non-coding predictions are
#' encoded as 1, and coding predictions as 0.
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
#' @return A list containing `seqIDs` (a character vector of all unique
#'   sequence IDs found across all provided files) and `tools` (a list where
#'   each element is a binary vector of predictions for a given tool).
#' @export
#' @importFrom utils read.table write.table read.csv2
#' @importFrom stats setNames
#' 
#' @examples
#' # --- 1. Create temporary output files for a reproducible example ---
#' cpc2File <- tempfile()
#' plekFile <- tempfile()
#'
#' # CPC2 format: ID ... classification
#' write.table(
#'   data.frame(V1="ID1", V8="noncoding"),
#'   cpc2File, col.names=FALSE, row.names=FALSE, sep="\t"
#' )
#' write.table(
#'   data.frame(V1="ID2", V8="coding"),
#'   cpc2File, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE
#' )
#'
#' # PLEK format: Classification Score >ID
#' write.table(
#'   data.frame(V1="Non-coding", V2="0.9", V3=">ID2"),
#'   plekFile, col.names=FALSE, row.names=FALSE, sep=" "
#' )
#'
#' # --- 2. Run the function with the temporary files ---
#' codpotResults <- aggregateCodPot(CPC2_outfile = cpc2File, PLEK_outfile = plekFile)
#'
#' print("Sequence IDs:")
#' print(codpotResults$seqIDs)
#' print("Tool predictions:")
#' print(codpotResults$tools)
#'
#' # --- 3. Clean up the temporary files ---
#' unlink(c(cpc2File, plekFile))
aggregateCodPot <- function(CPC2_outfile = NULL, PLEK_outfile = NULL,
                            FEELnc_outfile = NULL, CPAT_outfile = NULL,
                            CPAT_cutoff = 0.364, CNCI_outfile = NULL,
                            LncFinder_outfile = NULL,
                            lncRNA_Mdeep_outfile = NULL) {
    
    codpotList <- list(tools = list())
    
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
    
    codpotList$tools$CPC2 <- readFileIfExists(CPC2_outfile, readCpc2)
    codpotList$tools$PLEK <- readFileIfExists(PLEK_outfile, readPlek)
    codpotList$tools$FEELnc <- readFileIfExists(FEELnc_outfile, readFeelnc)
    codpotList$tools$CPAT <- readFileIfExists(CPAT_outfile, readCpat, CPAT_cutoff)
    codpotList$tools$CNCI <- readFileIfExists(CNCI_outfile, readCnci)
    codpotList$tools$LncFinder <- readFileIfExists(LncFinder_outfile, readLncFinder)
    codpotList$tools$lncRNA_Mdeep <- readFileIfExists(lncRNA_Mdeep_outfile,
                                                      readLncRnaMdeep)
    
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
        as.integer(seqIds %in% names(toolIds))
    })
    
    predictionMatrix <- do.call(cbind, toolPredictions)
    colsToKeep <- colSums(predictionMatrix, na.rm = TRUE) > 0
    resultList$tools <- toolPredictions[colsToKeep]
    
    return(resultList)
}

# --- Internal Helper Functions to Read Tool Outputs ---

#' @noRd
readCpc2 <- function(cpc2Outfile) {
    cpc2Data <- utils::read.table(cpc2Outfile, stringsAsFactors = FALSE)
    noncodingIds <- cpc2Data[cpc2Data$V8 == "noncoding", 1]
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readPlek <- function(plekOutfile) {
    plekData <- utils::read.table(plekOutfile, header = FALSE,
                                  stringsAsFactors = FALSE)
    if (ncol(plekData) < 3) {
        warning("PLEK output file '", plekOutfile,
                "' has fewer than 3 columns. Cannot extract IDs.")
        return(setNames(character(0), character(0)))
    }
    plekData$V3 <- gsub(">", "", plekData$V3)
    noncodingIds <- plekData[plekData$V1 == "Non-coding", 3]
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readFeelnc <- function(feelncOutfile) {
    feelncData <- utils::read.table(feelncOutfile, header = TRUE,
                                    stringsAsFactors = FALSE)
    noncodingIds <- feelncData[feelncData$label == 0, "name"]
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readCpat <- function(cpatOutfile, cpatCutoff) {
    cpatData <- utils::read.table(cpatOutfile, header = TRUE, row.names = 1,
                                  stringsAsFactors = FALSE)
    noncodingIds <- rownames(cpatData[cpatData$coding_prob < cpatCutoff, ,
                                      drop = FALSE])
    noncodingIds <- unique(noncodingIds)
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readCnci <- function(cnciOutfile) {
    cnciData <- utils::read.table(cnciOutfile, header = TRUE, sep = "\t",
                                  stringsAsFactors = FALSE)
    noncodingIds <- cnciData[cnciData$index == "noncoding", 1]
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readLncFinder <- function(lncFinderOutfile) {
    lncFinderData <- utils::read.csv2(lncFinderOutfile, header = TRUE,
                                      row.names = 1, stringsAsFactors = FALSE)
    noncodingIds <- rownames(lncFinderData[lncFinderData$Pred == "NonCoding", ,
                                           drop = FALSE])
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}

#' @noRd
readLncRnaMdeep <- function(lncRnaMdeepOutfile) {
    lncRnaMdeepData <- utils::read.table(lncRnaMdeepOutfile,
                                         stringsAsFactors = FALSE)
    noncodingIds <- lncRnaMdeepData[lncRnaMdeepData$V2 == "noncoding", 1]
    noncodingIds <- unique(as.character(noncodingIds))
    setNames(rep("noncoding", length(noncodingIds)), noncodingIds)
}


#' Create Venn Diagram from Coding Potential Results
#'
#' Generates a Venn diagram of noncoding predictions from multiple tools, based
#' on the output of `aggregateCodPot`. This function requires the 'venn' package,
#' which should be declared in the `Suggests` field of the DESCRIPTION file.
#'
#' @param codPot A list object generated by `aggregateCodPot()`.
#' @param selection An optional numeric vector indicating which tools to include
#'   (1 for include, 0 for exclude). If `NULL`, all available tools are used.
#' @param vennColors A vector of colors for the Venn diagram segments.
#'
#' @return Invisibly returns the result of `venn::venn`. The primary effect is
#'   plotting the diagram to the active graphics device.
#' @export
#'
#' @examples
#' # --- 1. Create a mock object that mimics the output of aggregateCodPot ---
#' mockCodPot <- list(
#'   seqIDs = c("ID1", "ID2", "ID3", "ID4"),
#'   tools = list(
#'     CPC2 = c(1, 1, 0, 0), # Noncoding: ID1, ID2
#'     PLEK = c(0, 1, 1, 0), # Noncoding: ID2, ID3
#'     CPAT = c(0, 0, 1, 1)  # Noncoding: ID3, ID4
#'   )
#' )
#'
#' # --- 2. Generate Venn diagrams (requires the 'venn' package) ---
#' if (requireNamespace("venn", quietly = TRUE)) {
#'   # To prevent plots from showing up during automated checks, we can
#'   # temporarily redirect the graphics output.
#'   png(tempfile())
#'
#'   # Example with all tools
#'   plotVennCodPot(codPot = mockCodPot)
#'
#'   # Example with a selection of tools (CPC2 and CPAT)
#'   plotVennCodPot(codPot = mockCodPot, selection = c(1, 0, 1))
#'
#'   dev.off() # Close the graphics device
#' } else {
#'   message("Skipping venn diagram examples, 'venn' package not installed.")
#' }
plotVennCodPot <- function(codPot, selection = NULL,
                           vennColors = c("green", "red", "blue", "yellow",
                                          "magenta", "black", "orange", "white",
                                          "darkgreen")) {
    
    if (!requireNamespace("venn", quietly = TRUE)) {
        message("The 'venn' package is required for this function.\n",
                "Please install it using: BiocManager::install('venn')")
        return(invisible(NULL))
    }
    
    stopifnot(
        "Input 'codPot' must be a list" = is.list(codPot),
        "'codPot' must contain 'seqIDs' and 'tools'" =
            all(c("seqIDs", "tools") %in% names(codPot)),
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
        stop("Length of 'selection' (", length(selection),
                   ") does not match the number of tools (", numTools, ").")
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
    
    # This function directly plots, so no object is returned from the call itself
    venn::venn(
        listForVenn,
        zcolor = vennColors[seq_len(numSelected)],
        snames = selectedToolNames,
        ilabels = TRUE
    )
    
    return(invisible(NULL))
}