#' Get Available Expression Units from StringTie Files
#'
#' Reads the first provided StringTie output file to determine the available
#' column names, which typically represent different expression units (e.g., FPKM, TPM).
#'
#' @param filePaths A character vector of full paths to one or more StringTie
#'   output files (typically with a `.gtab` extension).
#'
#' @return A character vector of column names from the first StringTie file.
#'
#' @export
#' @importFrom utils read.table write.table
#' @examples
#' # --- 1. Create a temporary StringTie file for the example ---
#' tempFile <- tempfile(fileext = ".gtab")
#' sampleContent <- data.frame(
#'   Gene.ID = "GENE1", Gene.Name = "G1", FPKM = 10.5, TPM = 12.1
#' )
#' write.table(sampleContent, tempFile, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # --- 2. Run the function ---
#' availableUnits <- getStringtieUnits(filePaths = tempFile)
#' print(availableUnits)
#'
#' # --- 3. Clean up ---
#' unlink(tempFile)
#'
getStringtieUnits <- function(filePaths) {
  stopifnot(
    "filePaths must be a character vector" = is.character(filePaths),
    "At least one file path must be provided" = (length(filePaths) > 0),
    "The first file does not exist" = file.exists(filePaths[1])
  )
  tab <- utils::read.table(filePaths[1], sep = "\t", header = TRUE, check.names = FALSE)
  return(colnames(tab))
}


#' Extract Expression Values from StringTie Files
#'
#' Reads multiple StringTie output files, extracts expression values for a
#' specified unit (e.g., FPKM or TPM), and combines them into a single data.frame.
#'
#' @param filePaths A character vector of full paths to StringTie output files.
#' @param expUnit A character string specifying the expression unit to extract
#'   (must be a column name in the files). Defaults to "FPKM".
#'
#' @return A data.frame with gene metadata and expression values, where each
#'   column represents a sample derived from the input file names.
#'
#' @export
#' @examples
#' # --- 1. Create a temporary directory and multiple StringTie files ---
#' tempDir <- tempdir()
#' file1 <- file.path(tempDir, "sampleA.gtab")
#' file2 <- file.path(tempDir, "sampleB.gtab")
#'
#' sampleContent1 <- data.frame(
#'   Gene.ID="G1", Gene.Name="GN1", Reference="chr1", Strand="+", Start=1, End=100, FPKM=10, TPM=12
#' )
#' sampleContent2 <- data.frame(
#'   Gene.ID="G1", Gene.Name="GN1", Reference="chr1", Strand="+", Start=1, End=100, FPKM=20, TPM=25
#' )
#' write.table(sampleContent1, file1, sep = "\t", row.names = FALSE, quote = FALSE)
#' write.table(sampleContent2, file2, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # --- 2. Run the function ---
#' filePaths <- c(file1, file2)
#' fpkmData <- getStringtieExpression(filePaths, expUnit = "FPKM")
#' print(fpkmData)
#'
#' # --- 3. Clean up ---
#' unlink(filePaths)
#'
getStringtieExpression <- function(filePaths, expUnit = "FPKM") {
  stopifnot(
    "filePaths must be a character vector" = is.character(filePaths),
    "At least one file path must be provided" = (length(filePaths) > 0),
    "One or more files do not exist" = all(file.exists(filePaths))
  )

  # Read all files into a list of data frames
  allData <- lapply(filePaths, utils::read.table, sep = "\t", header = TRUE, check.names = FALSE)

  # Validate that the expression unit exists in the first file
  if (!expUnit %in% colnames(allData[[1]])) {
    stop("Expression unit '", expUnit, "' not found in columns of file: ", filePaths[1])
  }

  # Extract metadata from the first file
  metadataCols <- c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start", "End")
  metadata <- allData[[1]][, intersect(metadataCols, colnames(allData[[1]]))]

  # Extract expression values from each file
  sampleNames <- sub("\\.gtab$", "", basename(filePaths))
  exprList <- lapply(seq_along(allData), function(i) {
    df <- allData[[i]]
    exprData <- data.frame(
      Gene.ID = df$"Gene.ID",
      value = df[[expUnit]],
      stringsAsFactors = FALSE
    )
    colnames(exprData)[2] <- sampleNames[i]
    return(exprData)
  })

  # Merge all expression data frames
  mergedExpr <- Reduce(function(x, y) merge(x, y, by = "Gene.ID", all = TRUE), exprList)

  # Join with metadata
  finalTable <- merge(metadata, mergedExpr, by = "Gene.ID", all.y = TRUE)
  return(finalTable)
}


#' Compute Gene Lengths from a StringTie File
#'
#' Computes gene lengths based on the 'Start' and 'End' coordinates from a
#' single StringTie output file.
#'
#' @param filePaths A character vector of full paths to one or more StringTie
#'   output files. Only the first file in the vector will be used.
#'
#' @return A data.frame with two columns: `Gene.ID` and `Length`.
#'
#' @export
#' @examples
#' # --- 1. Create a temporary StringTie file ---
#' tempFile <- tempfile(fileext = ".gtab")
#' sampleContent <- data.frame(
#'   Gene.ID = c("G1", "G2"), Start = c(100, 500), End = c(250, 600)
#' )
#' write.table(sampleContent, tempFile, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # --- 2. Run the function ---
#' geneLengths <- getStringtieGeneLength(filePaths = tempFile)
#' print(geneLengths)
#'
#' # --- 3. Clean up ---
#' unlink(tempFile)
#'
getStringtieGeneLength <- function(filePaths) {
  stopifnot(
    "filePaths must be a character vector" = is.character(filePaths),
    "At least one file path must be provided" = (length(filePaths) > 0),
    "The first file does not exist" = file.exists(filePaths[1])
  )
  # Gene length should be consistent, so we only need to read one file.
  tab <- utils::read.table(filePaths[1], sep = "\t", header = TRUE, check.names = FALSE)

  stopifnot(
    "Required columns 'Start' and 'End' not found" = all(c("Start", "End") %in% names(tab))
  )

  tab$Length <- abs(tab$End - tab$Start) + 1

  return(tab[, c("Gene.ID", "Length")])
}
