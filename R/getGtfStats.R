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
#' @importFrom dplyr %>% group_by summarise
#' @noRd
calculateTranscriptLength <- function(gtfObject) {
    transcript_id <- width <- NULL
    
    strGTF <- as.data.frame(gtfObject)
    strGTF_ex <- strGTF[strGTF$type %in% "exon", ]
    
    trans_len <- strGTF_ex %>%
        group_by(transcript_id) %>%
        summarise(trans_length = sum(width))
    
    trans_len <- as.data.frame(trans_len)
    return(trans_len)
}

#' Count Exons per Transcript (Internal Function)
#'
#' @param gtfObject A GRanges object.
#' @return A data.frame with transcript_id and exon count.
#' @importFrom dplyr %>% group_by summarise
#' @noRd
calculateExonCount <- function(gtfObject) {
    transcript_id <- exon_number <- NULL
    
    strGTF <- as.data.frame(gtfObject)
    strGTF_ex <- strGTF[strGTF$type %in% "exon", ]
    
    exons_n <- strGTF_ex %>%
        group_by(transcript_id) %>%
        summarise(exons = max(as.numeric(as.character(exon_number))))
    
    exons_n <- as.data.frame(exons_n)
    return(exons_n)
}