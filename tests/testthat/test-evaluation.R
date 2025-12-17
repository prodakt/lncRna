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
    
    # 3. Test BestTool
    tools_to_test <- c("ToolA", "ToolB")
    
    stats <- BestTool(summarySet, tools = tools_to_test)
    expect_s3_class(stats, "data.frame")
    expect_true("Accuracy" %in% rownames(stats))
    expect_equal(ncol(stats), 2)
    
    # 4. Test Thresholds
    threshSummary <- evaluateToolsThresholds(summarySet, tools = tools_to_test)
    expect_true("atLeastN" %in% names(threshSummary))
    expect_true("atl1" %in% names(threshSummary$atLeastN))
    
    threshStats <- BestToolAtleast(threshSummary)
    expect_s3_class(threshStats, "data.frame")
    
    # 5. Test Combinations
    combSummary <- evaluateToolCombinations(summarySet, tools = tools_to_test)
    
    expect_true("ToolA+ToolB" %in% names(combSummary$toolCombinations))
    
    combStats <- BestToolCombination(combSummary, combinations = "ToolA+ToolB")
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