# tests/testthat/test-interactions.R

library(testthat)
library(lncRna)

test_that("findTransInteractions filters correctly", {
  set.seed(123)
  n_samples <- 10
  expr <- matrix(rnorm(4 * n_samples), nrow = 4,
                 dimnames = list(c("LNC1", "T1", "LNC2", "T2"),
                                 paste0("S", 1:n_samples)))

  expr["LNC1", ] <- expr["T1", ] + rnorm(n_samples, 0, 0.01) # Positive corr
  expr["LNC2", ] <- -expr["T2", ] + rnorm(n_samples, 0, 0.01) # Negative corr

  res <- suppressMessages(findTransInteractions(
    expr,
    lncRnaList = c("LNC1", "LNC2"),
    tarRnaList = c("T1", "T2"),
    rval = 0.9, pval = 0.05
  ))

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

test_that("findCisInteractions works correctly with mock data", {
  feelnc_file <- tempfile(fileext = ".txt")

  mock_data <- data.frame(
    isBest = c(1, 1, 0, 1),
    lcRNA_transcript = c("LNC_T1", "LNC_T2", "LNC_T3", "LNC_T4"),
    lncRNA_gene = c("LNC_G1", "LNC_G2", "LNC_G3", "LNC_G4"),
    partnerRNA_gene = c("TARGET_G1", "TARGET_G2", "TARGET_G3", "TARGET_G4"),
    distance = c(5000, 150000, 1000, 9000),
    stringsAsFactors = FALSE
  )
  utils::write.table(mock_data, feelnc_file, sep = "\t", row.names = FALSE, col.names = TRUE)

  # Test 1: Domyślne filtrowanie (isBest == 1, maxDist <= 100000)
  res1 <- findCisInteractions(feelnc_file)

  expect_s3_class(res1, "data.frame")
  expect_equal(nrow(res1), 2)
  expect_true(all(res1$isBest == 1))
  expect_true(all(res1$distance <= 100000))

  expect_true("lncRNAId" %in% colnames(res1))
  expect_true("targetRNAId" %in% colnames(res1))

  # Test 2: Filtrowanie listą lncRNA i mRNA
  res2 <- findCisInteractions(
    FEELncClassesFile = feelnc_file,
    lncRnas = c("LNC_T1"),
    mRnas = c("TARGET_G1")
  )
  expect_equal(nrow(res2), 1)
  expect_equal(res2$lncRNAId, "LNC_G1")
  expect_equal(res2$targetRNAId, "TARGET_G1")

  unlink(feelnc_file)
})

test_that("findCisInteractions throws error for non-existent file", {
  expect_error(
    findCisInteractions("non_existent_feelnc_file.txt"),
    "The specified 'FEELncClassesFile' does not exist"
  )
})
