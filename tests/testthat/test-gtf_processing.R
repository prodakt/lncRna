test_that("getGtfStats correctly calculates statistics", {
    # Mock GRanges
    mockGtf <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = c(100, 200, 400), width = c(50, 60, 100)),
        type = c("exon", "exon", "transcript"),
        transcript_id = c("tx1", "tx1", "tx2"),
        exon_number = c("1", "2", NA)
    )
    
    stats <- getGtfStats(mockGtf)
    
    expect_s3_class(stats, "data.frame")
    expect_equal(nrow(stats), 1) # Only tx1 has exons
    expect_equal(stats$transcript_id, "tx1")
    expect_equal(stats$exons, 2)
    expect_equal(stats$trans_length, 110) # 50 + 60
})

test_that("getBiotypes extracts correct columns", {
    mockGtf <- GenomicRanges::GRanges(
        seqnames = "chr1", ranges = IRanges::IRanges(1, 10),
        gene_id = "G1", gene_biotype = "protein_coding",
        transcript_id = "T1", transcript_biotype = "lncRNA"
    )
    
    # Test gene level
    res_gene <- getBiotypes(mockGtf, level = "gene")
    expect_s4_class(res_gene, "DataFrame")
    expect_true("gene_id" %in% names(res_gene))
    expect_equal(res_gene$gene_biotype[1], "protein_coding")
    
    # Test transcript level
    res_trans <- getBiotypes(mockGtf, level = "transcript")
    expect_true("transcript_id" %in% names(res_trans))
    expect_equal(res_trans$transcript_biotype[1], "lncRNA")
})