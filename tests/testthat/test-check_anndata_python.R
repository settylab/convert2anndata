library(testthat)
library(convert2anndata)

test_that("check_anndata_python returns TRUE invisibly when env is fully working", {
  skip_if_not_installed("reticulate")
  if (!reticulate::py_module_available("anndata")) {
    skip("anndata module not available; the success path can't be exercised.")
  }
  withr::with_envvar(c(PYTHONPATH = ""), {
    expect_true(check_anndata_python())
  })
})

test_that("check_anndata_python action='silent' returns FALSE when anndata is missing", {
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py", version = "3.10"),
    py_module_available = function(name) FALSE,
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = ""), {
    expect_false(check_anndata_python(action = "silent"))
  })
})

test_that("check_anndata_python action='warn' produces a friendly warning", {
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py", version = "3.10"),
    py_module_available = function(name) FALSE,
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = ""), {
    expect_warning(
      res <- check_anndata_python(action = "warn"),
      regexp = "anndata.*not importable"
    )
    expect_false(res)
  })
})

test_that("check_anndata_python action='stop' raises with an actionable message", {
  local_mocked_bindings(
    py_available = function(...) FALSE,
    py_config    = function(...) list(python = "/fake/py", version = "3.10"),
    .package = "reticulate"
  )
  expect_error(
    check_anndata_python(),
    regexp = "No Python interpreter available"
  )
})

test_that("PYTHONPATH pointing at a different Python version is detected and reported", {
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py3.8/bin/python",
                                              version = "3.8.17"),
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = "/app/Python/3.11/site-packages"), {
    expect_error(
      check_anndata_python(),
      regexp = "PYTHONPATH appears to point at a different Python"
    )
  })
})

test_that("PYTHONPATH that matches the active Python version is allowed", {
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py3.10/bin/python",
                                              version = "3.10.12"),
    py_module_available = function(name) TRUE,
    py_run_string       = function(...) invisible(NULL),
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = "/some/path/python3.10/site-packages"), {
    expect_true(check_anndata_python())
  })
})

test_that("numpy import error gets a tailored hint", {
  fake_run <- function(code, ...) {
    if (grepl("import numpy", code, fixed = TRUE)) {
      stop("ImportError: Error importing numpy: you should not try to import numpy from its source directory")
    }
    invisible(NULL)
  }
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py", version = "3.10"),
    py_module_available = function(name) TRUE,
    py_run_string       = fake_run,
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = ""), {
    expect_error(
      check_anndata_python(),
      regexp = "PYTHONPATH is likely polluted"
    )
  })
})

test_that("anndata smoke-probe failure surfaces a tailored message", {
  fake_run <- function(code, ...) {
    if (grepl("anndata.AnnData", code, fixed = TRUE)) {
      stop("ValueError: numpy.dtype size changed (ABI break)")
    }
    invisible(NULL)
  }
  local_mocked_bindings(
    py_available        = function(...) TRUE,
    py_config           = function(...) list(python = "/fake/py", version = "3.10"),
    py_module_available = function(name) TRUE,
    py_run_string       = fake_run,
    .package = "reticulate"
  )
  withr::with_envvar(c(PYTHONPATH = ""), {
    expect_error(
      check_anndata_python(),
      regexp = "smoke probe failed"
    )
  })
})

test_that("diagnose_anndata_python prints a diagnostic block and returns logical", {
  out <- capture.output(res <- diagnose_anndata_python())
  expect_true(any(grepl("convert2anndata Python environment", out)))
  expect_true(any(grepl("python:", out)))
  expect_true(any(grepl("RETICULATE_PYTHON", out)))
  expect_type(res, "logical")
})

test_that("verbose=TRUE on success prints the env diagnostic block", {
  skip_if_not_installed("reticulate")
  if (!reticulate::py_module_available("anndata")) skip("anndata not available")
  withr::with_envvar(c(PYTHONPATH = ""), {
    out <- capture.output(check_anndata_python(verbose = TRUE))
    expect_true(any(grepl("convert2anndata Python environment", out)))
  })
})
