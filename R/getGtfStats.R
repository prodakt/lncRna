#' Extract Transcript Statistics from a GTF object
#'
#' This function takes a GRanges object (imported from a GTF file) and
#' calculates the number of exons and total transcript length for each
#' transcript.
#'
#' @param gtfObject A \code{GRanges} object, typically imported from a GTF
#'   file using \code{rtracklayer::import()}.
#'
#' @return A \code{data.frame} with columns: "transcript_id", "exons", and
#'   "trans_length".
#'
#' @export
#'
#' @examples
#' # Create a sample GRanges object to mimic a GTF import
#' sample_gtf <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(
#'     start = c(100, 300, 800, 950),
#'     end = c(200, 400, 900, 1050)
#'   ),
#'   strand = "+",
#'   type = "exon",
#'   transcript_id = c("T1", "T1", "T2", "T2"),
#'   exon_number = c("1", "2", "1", "20") # Example with exon_number > 9
#' )
#'
#' # Calculate statistics
#' transcript_stats <- getGtfStats(sample_gtf)
#' print(transcript_stats)
getGtfStats <- function(gtfObject) {
    stats <- merge(
        calculateExonCount(gtfObject),
        calculateTranscriptLength(gtfObject),
        by = "transcript_id",
        all = TRUE
    )
    return(stats)
}

#' Calculate Total Transcript Length (Internal Function)
#'
#' @param gtfObject A GRanges object.
#' @return A data.frame with transcript_id and trans_length.
#' @noRd
calculateTranscriptLength <- function(gtfObject) {
  strGTF <- as.data.frame(gtfObject)
  strGTF_ex <- strGTF[strGTF$type == "exon", ]

  if (nrow(strGTF_ex) == 0) {
    return(data.frame(transcript_id = character(),
                      trans_length = numeric()))
  }
  trans_len <- stats::aggregate(width ~ transcript_id,
                                data = strGTF_ex,
                                FUN = sum)

  colnames(trans_len)[2] <- "trans_length"

  return(trans_len)
}

#' Calculate Exon Count per Transcript
#'
#' @param gtfObject A GRanges object.
#' @return A data.frame with transcript_id and exons count.
#' @noRd
calculateExonCount <- function(gtfObject) {
  strGTF <- as.data.frame(gtfObject)
  strGTF_ex <- strGTF[strGTF$type == "exon", ]

  if (nrow(strGTF_ex) == 0) return(data.frame())
  exons_n <- stats::aggregate(type ~ transcript_id,
                              data = strGTF_ex,
                              FUN = length)

  colnames(exons_n)[2] <- "exons"
  return(exons_n)
}
