#' Get Gene or Transcript Biotypes from a GRanges Object
#'
#' This function extracts unique biotypes for either genes or transcripts from a
#' `GRanges` object, which is the standard representation for genomic features
#' in Bioconductor (e.g., after importing a GTF file with `rtracklayer::import`).
#'
#' @param gtfObject A `GRanges` object. The object's metadata columns must
#'   contain `gene_id` and `gene_biotype` for gene-level extraction, or
#'   `transcript_id` and `transcript_biotype` for transcript-level extraction.
#' @param level A character string specifying whether to extract information at
#'   the "gene" or "transcript" level. Defaults to "transcript".
#'
#' @return A `data.frame` with two columns containing unique identifier-biotype
#'   pairs (e.g., `transcript_id` and `transcript_biotype`), with rows
#'   containing `NA` identifiers removed.
#'
#' @keywords GTF biotype GRanges rtracklayer
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges GRanges
#' @export
#' @examples
#' # --- 1. Create a sample GRanges object for a reproducible example ---
#' # This object mimics the structure of a GTF file imported with rtracklayer.
#' sampleGtf <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = c(100, 200, 300), width = 50),
#'   strand = "+",
#'   gene_id = c("gene1", "gene1", NA),
#'   gene_biotype = c("protein_coding", "protein_coding", "lincRNA"),
#'   transcript_id = c("tx1", "tx2", "tx3"),
#'   transcript_biotype = c("protein_coding", "retained_intron", "lincRNA")
#' )
#'
#' # --- 2. Run the function ---
#'
#' # Example 2a: Extract transcript biotypes (default level)
#' transcriptBiotypes <- getRefBiotypes(gtfObject = sampleGtf)
#' print("Transcript Biotypes:")
#' print(transcriptBiotypes)
#'
#' # Example 2b: Extract gene biotypes
#' geneBiotypes <- getRefBiotypes(gtfObject = sampleGtf, level = "gene")
#' print("Gene Biotypes:")
#' print(geneBiotypes)
#'
getRefBiotypes <- function(gtfObject, level = "transcript") {

  # --- 1. Input Validation ---
  stopifnot(
    "'gtfObject' must be a GRanges object" = inherits(gtfObject, "GRanges"),
    "'level' must be either 'gene' or 'transcript'" = (level %in% c("gene", "transcript"))
  )

  # Convert only the metadata columns to a data.frame for efficiency
  metadataDf <- as.data.frame(S4Vectors::mcols(gtfObject))

  # --- 2. Extract Biotypes based on Level ---
  if (level == "gene") {
    requiredCols <- c("gene_id", "gene_biotype")
    if (!all(requiredCols %in% names(metadataDf))) {
      stop("Missing required metadata columns for level='gene': 'gene_id', 'gene_biotype'")
    }
    result <- unique(metadataDf[, requiredCols])
    result <- result[!is.na(result$gene_id), ]
  } else { # level == "transcript"
    requiredCols <- c("transcript_id", "transcript_biotype")
    if (!all(requiredCols %in% names(metadataDf))) {
      stop("Missing required metadata columns for level='transcript': 'transcript_id', 'transcript_biotype'")
    }
    result <- unique(metadataDf[, requiredCols])
    result <- result[!is.na(result$transcript_id), ]
  }

  return(result)
}
