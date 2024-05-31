library(testthat)
library(withr)

# Define the shell script path
shell_script <- "tests/testthat/run_cli_tests.sh"

# Helper function to run shell commands
run_shell_command <- function(cmd) {
  result <- system(cmd, intern = TRUE)
  return(result)
}

test_that("cli_convert works with SingleCellExperiment input", {
  result <- run_shell_command(paste("bash", shell_script, "SingleCellExperiment input"))
  expect_true(any(grepl("Test passed: SingleCellExperiment input", result)))
})

test_that("cli_convert works with Seurat input", {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    result <- run_shell_command(paste("bash", shell_script, "Seurat input"))
    expect_true(any(grepl("Test passed: Seurat input", result)))
  }
})

test_that("cli_convert uses default assay", {
  result <- run_shell_command(paste("bash", shell_script, "Default assay"))
  expect_true(any(grepl("Test passed: Default assay", result)))
})

test_that("cli_convert handles altExp flag", {
  result <- run_shell_command(paste("bash", shell_script, "altExp flag"))
  expect_true(any(grepl("Test passed: altExp flag", result)))
})

test_that("cli_convert stops without input", {
  result <- run_shell_command(paste("bash", shell_script, "No input"))
  expect_true(any(grepl("Test passed: No input", result)))
})
