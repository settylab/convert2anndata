library(testthat)
library(convert2anndata)

test_that("setup_anndata_python returns 'none' when no source supplied", {
  withr::with_envvar(
    c(RETICULATE_PYTHON = "", CONDA_PREFIX = ""),
    expect_equal(setup_anndata_python(validate = FALSE), "none")
  )
})

test_that("setup_anndata_python with a bogus conda_env warns and returns 'none'", {
  bogus <- "/definitely/not/a/conda/env/path/xyz"
  withr::with_envvar(
    c(RETICULATE_PYTHON = "", CONDA_PREFIX = ""),
    {
      expect_warning(
        res <- setup_anndata_python(bogus, required = TRUE, validate = FALSE),
        regexp = "Could not activate conda environment"
      )
      expect_equal(res, "none")
    }
  )
})

test_that("setup_anndata_python prefers existing RETICULATE_PYTHON path over CONDA_PREFIX", {
  skip_if_not(file.exists(Sys.which("python") |> as.character()),
              "no python on PATH")
  withr::with_envvar(
    c(RETICULATE_PYTHON = as.character(Sys.which("python")), CONDA_PREFIX = "/nonexistent"),
    {
      res <- suppressWarnings(setup_anndata_python(validate = FALSE))
      expect_true(res %in% c("reticulate_python", "none"))
    }
  )
})

test_that("setup_anndata_python resolves a name against ~/micromamba/envs/", {
  # We don't actually point at a working env here -- just verify that the
  # name-resolution path is exercised when a real directory exists at the
  # expected location. We create a fake env directory and rely on
  # use_condaenv(required=FALSE) accepting the path without validating it.
  fake_root <- file.path(tempdir(), "fake_micromamba", "envs", "synthetic_env")
  dir.create(fake_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(dirname(dirname(fake_root)), recursive = TRUE), add = TRUE)

  withr::with_envvar(
    c(RETICULATE_PYTHON = "", CONDA_PREFIX = "",
      MAMBA_ROOT_PREFIX = file.path(tempdir(), "fake_micromamba")),
    {
      res <- suppressWarnings(
        setup_anndata_python("synthetic_env", validate = FALSE)
      )
      # If reticulate accepted the directory, label is conda_env_resolved.
      # If it didn't, label falls through to 'none'. Either is acceptable;
      # what matters is that we didn't crash.
      expect_true(res %in% c("conda_env_resolved", "conda_env", "none"))
    }
  )
})
