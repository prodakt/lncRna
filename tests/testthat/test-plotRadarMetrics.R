library(testthat)
library(lncRna)

test_that("plotRadarMetrics runs without errors", {
  mockCmList <- list(
    `MethodA` = list(metrics = c(Accuracy=0.9, Sensitivity=0.8, Specificity=0.95, Precision=0.85, Recall=0.8)),
    `MethodB` = list(metrics = c(Accuracy=0.8, Sensitivity=0.9, Specificity=0.7, Precision=0.75, Recall=0.9))
  )

  pdf(file = NULL)

  expect_invisible(plotRadarMetrics(cmList = mockCmList, methods = c("MethodA", "MethodB")))

  expect_invisible(plotRadarMetrics(cmList = mockCmList, methods = c("MethodA", "MethodB"),
                                    layout = "multiple", displayArea = TRUE))

  dev.off()
})

test_that("plotRadarMetrics saves data to csv correctly", {
  mockCmList <- list(
    `MethodA` = list(metrics = c(Accuracy=0.9, Sensitivity=0.8, Specificity=0.95, Precision=0.85, Recall=0.8))
  )

  temp_base <- tempfile()

  pdf(file = NULL)
  expect_message(
    plotRadarMetrics(cmList = mockCmList, methods = "MethodA",
                     saveData = TRUE, fileName = temp_base),
    "Radar data frame saved to"
  )
  dev.off()

  expected_file <- paste0(temp_base, ".csv")
  expect_true(file.exists(expected_file))

  unlink(expected_file)
})
