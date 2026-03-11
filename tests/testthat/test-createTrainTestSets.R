library(testthat)
library(lncRna)

test_that("createTrainTestSets splits character vectors correctly", {
  seqs <- paste0("seq", 1:10)

  set.seed(42)
  res <- createTrainTestSets(seqs, percentTrain = 0.8, prefix = "test_set")

  expect_type(res, "list")
  expect_named(res, c("test_set.train", "test_set.test"))

  expect_equal(length(res$test_set.train), 8)
  expect_equal(length(res$test_set.test), 2)

  expect_setequal(c(res$test_set.train, res$test_set.test), seqs)
})

test_that("createTrainTestSets handles named lists", {
  seqs_list <- list(a = 1, b = 2, c = 3, d = 4)
  set.seed(42)
  res <- createTrainTestSets(seqs_list, percentTrain = 0.5)

  expect_equal(length(res$set.train), 2)
  expect_true(all(res$set.train %in% c("a", "b", "c", "d")))
})

test_that("createTrainTestSets handles edge cases and errors", {
  expect_warning(
    res <- createTrainTestSets(character(0)),
    "No sequence names found"
  )
  expect_equal(length(res$set.train), 0)

  expect_error(
    createTrainTestSets(1:10),
    "character vector or an object with names"
  )

  expect_error(
    createTrainTestSets(paste0("seq", 1:10), percentTrain = 1.5),
    "number between 0 and 1"
  )
})
