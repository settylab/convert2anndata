#!/bin/bash

# Function to run a test and check if the output file is created
run_test() {
  cmd=$1
  output_file=$2
  expected_msg=$3

  echo "Running: $cmd"
  $cmd

  if [ -f $output_file ]; then
    echo "Test passed: $expected_msg"
    rm $output_file # Clean up
  else
    echo "Test failed: $expected_msg"
    exit 1
  fi
}

# Paths to the mock data files
mock_sce="tests/testdata/mock_sce.rds"
mock_seurat="tests/testdata/mock_seurat.rds"

# Ensure the mock data files exist
if [ ! -f $mock_sce ]; then
  echo "Mock SCE file not found: $mock_sce"
  exit 1
fi

if [ ! -f $mock_seurat ]; then
  echo "Mock Seurat file not found: $mock_seurat"
  exit 1
fi

case $1 in
  "SingleCellExperiment input")
    run_test "Rscript -e 'convert2anndata::cli_convert()' -i $mock_sce -o tests/testdata/output_sce.h5ad" "tests/testdata/output_sce.h5ad" "SingleCellExperiment input"
    ;;
  "Seurat input")
    run_test "Rscript -e 'convert2anndata::cli_convert()' -i $mock_seurat -o tests/testdata/output_seurat.h5ad" "tests/testdata/output_seurat.h5ad" "Seurat input"
    ;;
  "Default assay")
    run_test "Rscript -e 'convert2anndata::cli_convert()' -i $mock_sce -o tests/testdata/output_default_assay.h5ad -a counts" "tests/testdata/output_default_assay.h5ad" "Default assay"
    ;;
  "altExp flag")
    run_test "Rscript -e 'convert2anndata::cli_convert()' -i $mock_sce -o tests/testdata/output_altExp.h5ad -d" "tests/testdata/output_altExp.h5ad" "altExp flag"
    ;;
  "No input")
    run_test "Rscript -e 'convert2anndata::cli_convert()' -o tests/testdata/output_no_input.h5ad" "tests/testdata/output_no_input.h5ad" "No input"
    ;;
  *)
    echo "Invalid test case"
    exit 1
    ;;
esac
