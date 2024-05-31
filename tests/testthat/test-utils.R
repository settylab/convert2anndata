library(testthat)
library(convert2anndata)

test_that("timestamped_cat works correctly", {
  expect_output(timestamped_cat("Hello, world!"), "\\[\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}\\] Hello, world!")
})
