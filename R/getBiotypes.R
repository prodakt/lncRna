#' Extract Gene or Transcript Biotypes from a GRanges Object
#'
#' This function extracts unique biotypes for genes or transcripts from a
#' `GRanges` object, typically imported from a reference GTF/GFF file.
#'
#' @param refGtf A `GRanges` object imported from a reference annotation file.
#'   It must contain metadata columns for biotypes (e.g., `gene_biotype`,
#'   `transcript_biotype`) and corresponding IDs (`gene_id`, `transcript_id`).
#' @param level A character string specifying whether to extract information for
#'   `"gene"` or `"transcript"` features (default: `"transcript"`).
#'
#' @return A `DataFrame` (from the `S4Vectors` package) with two columns:
#'   the identifier (`gene_id` or `transcript_id`) and the corresponding
#'   biotype (`gene_biotype` or `transcript_biotype`).
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#' # --- 1. Create a sample GRanges object mimicking a reference GTF ---
#' sampleRefGtf <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1:6, width = 100),
#'   gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"),
#'   gene_biotype = c("protein_coding", "protein_coding", "lncRNA",
#'                    "lncRNA", "pseudogene", "pseudogene"),
#'   transcript_id = c("T1.1", "T1.2", "T2.1", "T2.1", "T3.1", "T3.1"),
#'   transcript_biotype = c("protein_coding", "protein_coding_variant",
#'                          "lncRNA", "lncRNA", "pseudogene", "pseudogene")
#' )
#'
#' # --- 2. Extract biotypes at the transcript level ---
#' transcript_biotypes <- getBiotypes(refGtf = sampleRefGtf, level = "transcript")
#' print(transcript_biotypes)
#'
#' # --- 3. Extract biotypes at the gene level ---
#' gene_biotypes <- getBiotypes(refGtf = sampleRefGtf, level = "gene")
#' print(gene_biotypes)
#'
getBiotypes <- function(refGtf, level = "transcript") {
    level <- match.arg(level, choices = c("transcript", "gene"))
    
    # Access the metadata columns directly from the GRanges object
    mcols_data <- S4Vectors::mcols(refGtf)
    
    if (level == "gene") {
        # Check for required columns
        if (!all(c("gene_id", "gene_biotype") %in% names(mcols_data))) {
            stop("Metadata columns 'gene_id' and 'gene_biotype' not found in refGtf.")
        }
        
        # Create a DataFrame with the two columns of interest
        biotypes_df <- S4Vectors::DataFrame(
            gene_id = mcols_data$gene_id,
            gene_biotype = mcols_data$gene_biotype
        )
        # Remove rows with NA IDs and then find unique combinations
        biotypes_df <- biotypes_df[!is.na(biotypes_df$gene_id), ]
        unique_biotypes <- unique(biotypes_df)
        
    } else { # level == "transcript"
        # Check for required columns
        if (!all(c("transcript_id", "transcript_biotype") %in% names(mcols_data))) {
            stop("Metadata columns 'transcript_id' and 'transcript_biotype' not found in refGtf.")
        }
        
        biotypes_df <- S4Vectors::DataFrame(
            transcript_id = mcols_data$transcript_id,
            transcript_biotype = mcols_data$transcript_biotype
        )
        biotypes_df <- biotypes_df[!is.na(biotypes_df$transcript_id), ]
        unique_biotypes <- unique(biotypes_df)
    }
    
    return(unique_biotypes)
}