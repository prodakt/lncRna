test_that("Evaluation pipeline works end-to-end", {
    # 1. Mock Data
    mockCodPotList <- list(
        seqIDs = c("tx1", "tx2", "tx3", "tx4"),
        tools = list(
            ToolA = c(0, 1, 1, 0),
            ToolB = c(0, 1, 0, 0)
        )
    )
    ncTest <- c("tx2", "tx3")
    cdsTest <- c("tx1", "tx4")

    # 2. Test prepareEvaluationSets
    summarySet <- prepareEvaluationSets(mockCodPotList, ncTest, cdsTest)
    expect_type(summarySet, "list")
    expect_equal(length(summarySet$seqIDs), 4)
    expect_equal(summarySet$type, c("cds", "nc", "nc", "cds"))

    # 3. Test bestTool
    tools_to_test <- c("ToolA", "ToolB")

    stats <- bestTool(summarySet, tools = tools_to_test)
    expect_s3_class(stats, "data.frame")
    expect_true("Accuracy" %in% rownames(stats))
    expect_equal(ncol(stats), 2)

    # 4. Test Thresholds
    threshSummary <- evaluateToolsThresholds(summarySet, tools = tools_to_test)
    expect_true("atLeastN" %in% names(threshSummary))
    expect_true("atl1" %in% names(threshSummary$atLeastN))

    threshStats <- bestToolAtleast(threshSummary)
    expect_s3_class(threshStats, "data.frame")

    # 5. Test Combinations
    combSummary <- evaluateToolCombinations(summarySet, tools = tools_to_test)

    expect_true("ToolA+ToolB" %in% names(combSummary$toolCombinations))

    combStats <- bestToolCombination(combSummary, combinations = "ToolA+ToolB")
    expect_s3_class(combStats, "data.frame")
})

test_that("calculateCM returns correct structure", {
    mockSummary <- list(
        isNC = c(1, 0, 1, 0),
        toolCombinations = list(Combo1 = c(1, 0, 0, 0))
    )

    cms <- calculateCM(mockSummary)
    expect_type(cms, "list")
    expect_true("Combo1" %in% names(cms))
    expect_true("table" %in% names(cms$Combo1))
    expect_true("metrics" %in% names(cms$Combo1))
    expect_true("positive" %in% names(cms$Combo1))
})

test_that("calculateCM handles filtering and printing correctly", {
  # 1.
  set.seed(123)
  mockSummary <- list(
    isNC = c(1, 1, 0, 0),
    toolCombinations = list(
      `Combo1` = c(1, 1, 0, 0),
      `Combo2` = c(0, 0, 1, 1)
    )
  )

  # 2.
  res_all <- calculateCM(mockSummary)
  expect_length(res_all, 2)
  expect_named(res_all, c("Combo1", "Combo2"))
  expect_true("metrics" %in% names(res_all$Combo1))

  # 3.
  res_high <- calculateCM(mockSummary, threshold = 0.8, returnOnlyHigh = TRUE)
  expect_length(res_high, 1)
  expect_named(res_high, "Combo1")

  # 4.
  expect_message(
    calculateCM(mockSummary, printMetricThresholds = TRUE),
    "Confusion Matrix Calculation & Filtering"
  )

  # 5.
  mockMetrics <- data.frame(
    Combo1 = c(Accuracy = 1.0, Sensitivity = 1.0),
    Combo2 = c(Accuracy = 0.0, Sensitivity = 0.0),
    row.names = c("Accuracy", "Sensitivity")
  )

  res_metrics <- calculateCM(
    mockSummary,
    metricsData = mockMetrics,
    printMetricThresholds = TRUE,
    metricToExtract = "Sensitivity"
  )
  expect_length(res_metrics, 2)
})

test_that("calculateCM handles edge cases", {
  expect_warning(
    calculateCM(list(isNC = c(0,1), toolCombinations = list())),
    "is empty"
  )

  expect_warning(
    calculateCM(list(isNC = c(1,1), toolCombinations = list(A=c(1,1)))),
    "must contain both classes"
  )
})
