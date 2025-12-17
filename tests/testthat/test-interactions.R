test_that("findTransInteractions filters correctly", {
    set.seed(123)
    n_samples <- 10
    expr <- matrix(rnorm(4 * n_samples), nrow = 4,
                   dimnames = list(c("LNC1", "T1", "LNC2", "T2"), 
                                   paste0("S", 1:n_samples)))
    
    expr["LNC1", ] <- expr["T1", ] + rnorm(n_samples, 0, 0.01) # Positive corr
    expr["LNC2", ] <- -expr["T2", ] + rnorm(n_samples, 0, 0.01) # Negative corr
    
    res <- findTransInteractions(
        expr, 
        lncRnaList = c("LNC1", "LNC2"), 
        tarRnaList = c("T1", "T2"),
        rval = 0.9, pval = 0.05
    )
    
    expect_s3_class(res, "data.frame")
    found_pos <- any((res$lncRNAId == "LNC1" & res$targetRNAId == "T1") | 
                         (res$lncRNAId == "T1" & res$targetRNAId == "LNC1"))
    expect_true(found_pos)
})

test_that("annotateInteractions merges correctly", {
    gostRes <- data.frame(
        term_name = "Stress", 
        intersection = "T1,T2",
        stringsAsFactors = FALSE
    )
    interactions <- data.frame(
        lncRNAId = c("L1", "L2"),
        targetRNAId = c("T1", "T3"), 
        stringsAsFactors = FALSE
    )
    
    res <- annotateInteractions(
        gostResult = gostRes, 
        interactionTable = interactions, 
        type = "cis",
        lncRnaCol = "lncRNAId",
        targetCol = "targetRNAId"
    )
    
    expect_s3_class(res, "data.frame")
    expect_true("L1" %in% res$lncRNAId)
    expect_equal(unique(res$type), "cis")
})