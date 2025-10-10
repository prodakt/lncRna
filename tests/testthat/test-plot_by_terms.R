library(testthat)

# Test 2: Multiple rows with color_selected
test_that("plot_by_terms handles multiple rows and color_selected", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("Polychrome")

  data <- data.frame(
    term_name = c("term1", "term2"),
    intersection = c("gene1", "gene2"),
    lncRNA.id = c("lnc1", "lnc2"),
    type = c("trans", "cis"),
    stringsAsFactors = FALSE
  )
  fig <- plot_by_terms(data, select_terms = "term1", color_selected = TRUE)
  expect_true(inherits(fig, "plotly"))
})

# Test 3: Error on missing columns
test_that("plot_by_terms errors with missing required columns", {
  data <- data.frame(term_name = "term1", lncRNA.id = "lnc1")  # Missing intersection, type
  expect_error(plot_by_terms(data), "Input 'data' must contain columns")
})

# Test 4: Error on invalid select_terms
test_that("plot_by_terms errors with invalid select_terms", {
  data <- data.frame(
    term_name = "term1",
    intersection = "gene1",
    lncRNA.id = "lnc1",
    type = "trans",
    stringsAsFactors = FALSE
  )
  expect_error(plot_by_terms(data, select_terms = "termX"), "Selected terms not present")
})

# Test 5: Error on color_selected without select_terms
test_that("plot_by_terms errors with color_selected but no select_terms", {
  data <- data.frame(
    term_name = "term1",
    intersection = "gene1",
    lncRNA.id = "lnc1",
    type = "trans",
    stringsAsFactors = FALSE
  )
  expect_error(plot_by_terms(data, color_selected = TRUE), "provide 'select_terms'")
})
