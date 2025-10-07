#' Extract Transcript Statistics from a GRanges Object
#'
#' This function calculates per-transcript statistics, such as the number of exons
#' and total transcript length (sum of exon lengths), from a `GRanges` object
#' typically created by importing a GTF/GFF file.
#'
#' @param gtfObject A `GRanges` object, for example, from `rtracklayer::import`.
#'   The object's metadata columns must contain `type`, `transcript_id`, and
#'   `exon_number` for the features of type "exon".
#'
#' @return A `data.frame` with columns: `transcript_id`, `transLength` (total
#'   length of all exons for the transcript), and `exons` (the highest exon
#'   number, representing the count of exons).
#'
#' @keywords GTF stringtie lncRNA GRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges GRanges width split
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' # --- 1. Create a sample GRanges object for a reproducible example ---
#' # This object mimics the structure of a GTF file imported with rtracklayer.
#' sampleGtf <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = c(100, 200, 400, 500), width = c(50, 60, 100, 20)),
#'   strand = "+",
#'   type = c("exon", "exon", "transcript", "exon"),
#'   transcript_id = c("tx1", "tx1", "tx2", "tx2"),
#'   exon_number = c("1", "2", NA, "1")
#' )
#'
#' # --- 2. Run the function ---
#' transcriptStats <- getGtfStats(gtfObject = sampleGtf)
#'
#' print("Transcript Statistics:")
#' print(transcriptStats)
#' # Expected output will show stats for tx1 (length=110, exons=2) and tx2 (length=20, exons=1).
#'
getGtfStats <- function(gtfObject) {

  # --- 1. Input Validation ---
  stopifnot(
    "'gtfObject' must be a GRanges object" = inherits(gtfObject, "GRanges")
  )

  # --- 2. Filter for Exons (Efficient Bioconductor way) ---
  # Work directly on the GRanges object, avoiding conversion to data.frame
  metadata <- S4Vectors::mcols(gtfObject)
  exons <- gtfObject[metadata$type == "exon"]

  if (length(exons) == 0) {
    warning("No features of type 'exon' found in the gtfObject.")
    return(data.frame(
      transcript_id = character(0),
      transLength = integer(0),
      exons = integer(0)
    ))
  }

  exonsMetadata <- S4Vectors::mcols(exons)
  requiredCols <- c("transcript_id", "exon_number")
  if (!all(requiredCols %in% names(exonsMetadata))) {
    stop("Metadata for exons is missing required columns: ",
         paste(setdiff(requiredCols, names(exonsMetadata)), collapse = ", "))
  }

  # --- 3. Calculate Stats using GRangesList (very fast) ---
  # Split the GRanges object into a list of GRanges, one for each transcript
  exonsByTx <- GenomicRanges::split(exons, exonsMetadata$transcript_id)

  # Calculate total length of exons per transcript (vectorized operation)
  transLength <- sum(GenomicRanges::width(exonsByTx))

  # Calculate number of exons per transcript (highest exon number)
  exonCount <- vapply(exonsByTx, function(tx) {
    # Ensure exon_number is numeric for max()
    nums <- as.numeric(S4Vectors::mcols(tx)$exon_number)
    if (all(is.na(nums))) return(NA_integer_)
    max(nums, na.rm = TRUE)
  }, integer(1))

  # --- 4. Combine Results into a Final data.frame ---
  statsDf <- data.frame(
    transcript_id = names(transLength),
    transLength = as.integer(transLength),
    exons = as.integer(exonCount),
    row.names = NULL
  )

  return(statsDf)
}
