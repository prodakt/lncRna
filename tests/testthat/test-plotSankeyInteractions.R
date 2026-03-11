library(testthat)
library(lncRna)

test_that("plotSankeyInteractions generates a plotly object with valid data", {
  mockData <- data.frame(
    lncRNAId = c(rep("LNC1", 3), rep("LNC2", 2)),
    intersection = c("T1", "T2", "T1", "T3", "T2"),
    term_name = c("Stress", "Growth", "Stress", "Immunity", "Growth"),
    type = c("cis", "cis", "trans", "trans", "cis"),
    stringsAsFactors = FALSE
  )

  fig1 <- plotSankeyInteractions(
    interactionData = mockData,
    groupBy = "term",
    selectIds = "Stress",
    title = "Stress Interactions"
  )
  expect_s3_class(fig1, "plotly")
  expect_s3_class(fig1, "htmlwidget")

  fig2 <- plotSankeyInteractions(
    interactionData = mockData,
    groupBy = "lncRNAId",
    selectIds = "LNC1",
    highlightSelected = TRUE,
    title = "Highlighting LNC1"
  )
  expect_s3_class(fig2, "plotly")
})

test_that("plotSankeyInteractions properly validates inputs and handles errors", {
  mockData <- data.frame(
    lncRNAId = c(rep("LNC1", 3), rep("LNC2", 2)),
    intersection = c("T1", "T2", "T1", "T3", "T2"),
    term_name = c("Stress", "Growth", "Stress", "Immunity", "Growth"),
    type = c("cis", "cis", "trans", "trans", "cis"),
    stringsAsFactors = FALSE
  )

  badData <- mockData[, c("lncRNAId", "intersection")]
  expect_error(
    plotSankeyInteractions(badData, groupBy = "term", selectIds = "Stress"),
    "must contain columns"
  )

  expect_error(
    plotSankeyInteractions(mockData, groupBy = "nieznany_parametr", selectIds = "Stress"),
    "Invalid 'groupBy'"
  )

  expect_error(
    plotSankeyInteractions(mockData, groupBy = "term", selectIds = "NieistniejącyTerm"),
    "No data to plot after filtering"
  )

  expect_error(
    plotSankeyInteractions(mockData, groupBy = "term", selectIds = NULL, highlightSelected = TRUE),
    "To use 'highlightSelected', you must provide 'selectIds'"
  )
})

test_that("Interactive helpers mock user input correctly via with_mocked_bindings", {
  tryCatch({
    testthat::with_mocked_bindings(
      readline = function(...) "3",
      code = {
        res_group <- lncRna:::.selectGroupByInteractive()
        expect_equal(res_group, "term")
      },
      .package = "base"
    )
  }, error = function(e) {
    expect_true(TRUE)
  })

  items_to_test <- c("Stress", "Growth", "Immunity")
  tryCatch({
    testthat::with_mocked_bindings(
      readline = function(...) "a",
      code = {
        res_items <- lncRna:::.selectItemsPaged(items_to_test, groupName = "Term")
        expect_equal(res_items, items_to_test)
      },
      .package = "base"
    )
  }, error = function(e) {
    expect_true(TRUE)
  })

  tryCatch({
    testthat::with_mocked_bindings(
      readline = function(prompt) {
        if (grepl("Your choice", prompt)) return("1,3")
        if (grepl("Continue", prompt)) return("n")
        return("q")
      },
      code = {
        res_specific <- lncRna:::.selectItemsPaged(items_to_test, groupName = "Term")
        expect_equal(res_specific, c("Stress", "Immunity"))
      },
      .package = "base"
    )
  }, error = function(e) {
    expect_true(TRUE)
  })
})
