library(testthat)
library(lncRna)

test_that("plotClockMetrics generates correct plot objects", {
  mockCmList <- list(
    `MethodA` = list(metrics = c(Accuracy=0.9, Sensitivity=0.8, Specificity=0.95, Precision=0.85, Recall=0.8)),
    `MethodB` = list(metrics = c(Accuracy=0.8, Sensitivity=0.9, Specificity=0.7, Precision=0.75, Recall=0.9))
  )

  p1 <- plotClockMetrics(
    cmList = mockCmList,
    methods = c("MethodA", "MethodB"),
    layout = "single"
  )
  expect_s3_class(p1, "ggplot")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p2 <- plotClockMetrics(
      cmList = mockCmList,
      methods = c("MethodA", "MethodB"),
      layout = "multiple"
    )
    expect_s3_class(p2, "patchwork")
  }

  expect_error(
    plotClockMetrics(cmList = mockCmList, methods = "MethodA", metrics = c("FakeMetric")),
    "Invalid metrics specified"
  )

  expect_error(
    plotClockMetrics(cmList = list()),
    "must be a non-empty, named list"
  )
})
