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
#'
getGtfStats <- function(gtfObject) {

  stopifnot(
    "'gtfObject' must be a GRanges object" = inherits(gtfObject, "GRanges")
  )

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

  exonsByTx <- GenomicRanges::split(exons, exonsMetadata$transcript_id)

  transLength <- sum(GenomicRanges::width(exonsByTx))

  exonCount <- vapply(exonsByTx, function(tx) {
    nums <- as.numeric(S4Vectors::mcols(tx)$exon_number)
    if (all(is.na(nums))) return(NA_integer_)
    as.integer(max(nums, na.rm = TRUE))
  }, integer(1))

  statsDf <- data.frame(
    transcript_id = names(transLength),
    transLength = as.integer(transLength),
    exons = as.integer(exonCount),
    row.names = NULL
  )

  return(statsDf)
}
