# tests/testthat/test-aggregateCodPot.R

library(testthat)
library(lncRna)

test_that("aggregateCodPot handles missing or empty inputs gracefully", {
  expect_warning(
    res <- aggregateCodPot(),
    "No valid tool output files were provided or found."
  )
  expect_type(res, "list")
  expect_length(res$seqIDs, 0)
  expect_length(res$tools, 0)
})

test_that("aggregateCodPot correctly parses mock tool outputs", {
  cpc2_file <- tempfile(fileext = ".txt")
  plek_file <- tempfile(fileext = ".txt")

  cat("tx1\t0\t0\t0\t0\t0\t0\tnoncoding\n",
      "tx2\t0\t0\t0\t0\t0\t0\tcoding\n",
      "tx3\t0\t0\t0\t0\t0\t0\tnoncoding\n",
      "tx4\t0\t0\t0\t0\t0\t0\tcoding\n",
      file = cpc2_file, sep = "")

  cat("Coding\t0.9\t>tx1\n",
      "Non-coding\t0.9\t>tx2\n",
      "Coding\t0.9\t>tx3\n",
      "Non-coding\t0.9\t>tx4\n",
      file = plek_file, sep = "")

  res <- aggregateCodPot(CPC2_outfile = cpc2_file, PLEK_outfile = plek_file)

  expect_type(res, "list")

  expect_true("CPC2" %in% names(res$tools))
  expect_true("PLEK" %in% names(res$tools))
  expect_true(all(c("tx1", "tx2", "tx3", "tx4") %in% res$seqIDs))
  cpc2_named <- setNames(res$tools$CPC2, res$seqIDs)
  expect_equal(as.numeric(cpc2_named[c("tx1", "tx2", "tx3", "tx4")]), c(1, 0, 1, 0))

  plek_named <- setNames(res$tools$PLEK, res$seqIDs)
  expect_equal(as.numeric(plek_named[c("tx1", "tx2", "tx3", "tx4")]), c(0, 1, 0, 1))

  unlink(c(cpc2_file, plek_file))
})

test_that("plotVennCodPot handles basic validations", {
  mockCodPot <- list(
    seqIDs = c("tx1", "tx2", "tx3", "tx4"),
    tools = list(
      CPC2 = c(1, 1, 0, 0),
      PLEK = c(0, 1, 1, 0)
    )
  )

  expect_error(plotVennCodPot(codPot = c(1, 2, 3)), "Input 'codPot' must be a list")
  expect_error(plotVennCodPot(codPot = mockCodPot, selection = c(TRUE)), "does not match the number of tools")

  if (requireNamespace("venn", quietly = TRUE)) {
    pdf(file = NULL)
    plotVennCodPot(codPot = mockCodPot)
    dev.off()
    expect_true(TRUE)
  }
})
