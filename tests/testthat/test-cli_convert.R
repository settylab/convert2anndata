library(testthat)
library(withr)

# Function to run shell commands and capture output
run_shell_command <- function(cmd) {
  cat("Running: ", cmd, "\n")
  result <- system(cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  return(result)
}

shell_script <- "./run_cli_tests.sh"
cat("\nCurrent directory:", getwd(), "\n")
cat("Bash script dir is: ", shell_script, "\n")

test_that("cli_convert works with SingleCellExperiment input", {
  result <- run_shell_command(paste("bash", shell_script, "SingleCellExperiment", "input"))
  cat("Result:\n", paste(result, collapse = "\n"), "\n")
  expect_true(any(grepl("Test passed: SingleCellExperiment input", result)))
})

if (requireNamespace("Seurat", quietly = TRUE)) {
  test_that("cli_convert works with Seurat input", {
    result <- run_shell_command(paste("bash", shell_script, "Seurat", "input"))
    cat("Result:\n", paste(result, collapse = "\n"), "\n")
    expect_true(any(grepl("Test passed: Seurat input", result)))
  })
}

test_that("cli_convert uses default assay", {
  result <- run_shell_command(paste("bash", shell_script, "Default", "assay"))
  cat("Result:\n", paste(result, collapse = "\n"), "\n")
  expect_true(any(grepl("Test passed: Default assay", result)))
})

test_that("cli_convert handles altExp flag", {
  result <- run_shell_command(paste("bash", shell_script, "altExp", "flag"))
  cat("Result:\n", paste(result, collapse = "\n"), "\n")
  expect_true(any(grepl("Test passed: altExp flag", result)))
})

test_that("cli_convert stops without input", {
  result <- run_shell_command(paste("bash", shell_script, "No", "input"))
  cat("Result:\n", paste(result, collapse = "\n"), "\n")
  expect_true(any(grepl("Test passed: No input", result)))
})
