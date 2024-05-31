#!/bin/bash

test_case=$1
test_input=$2

# Determine the directory of the current script
script_dir="$(dirname "$(realpath "$0")")"
mock_data_dir="$script_dir/../testdata"

case $test_case in
  "SingleCellExperiment")
    echo "Running SingleCellExperiment input test"
    Rscript -e "convert2anndata::cli_convert()" -i "$mock_data_dir/mock_sce.rds" -o "$mock_data_dir/output_sce.h5ad"
    if [ -f "$mock_data_dir/output_sce.h5ad" ]; then
      echo "Test passed: SingleCellExperiment input"
    else
      echo "Test failed: SingleCellExperiment input"
    fi
    ;;
  "Seurat")
    echo "Running Seurat input test"
    Rscript -e "convert2anndata::cli_convert()" -i "$mock_data_dir/mock_seurat.rds" -o "$mock_data_dir/output_seurat.h5ad"
    if [ -f "$mock_data_dir/output_seurat.h5ad" ]; then
      echo "Test passed: Seurat input"
    else
      echo "Test failed: Seurat input"
    fi
    ;;
  "Default")
    echo "Running default assay test"
    Rscript -e "convert2anndata::cli_convert()" -i "$mock_data_dir/mock_sce.rds" -o "$mock_data_dir/output_default_assay.h5ad" -a "counts"
    if [ -f "$mock_data_dir/output_default_assay.h5ad" ]; then
      echo "Test passed: Default assay"
    else
      echo "Test failed: Default assay"
    fi
    ;;
  "altExp")
    echo "Running altExp flag test"
    Rscript -e "convert2anndata::cli_convert()" -i "$mock_data_dir/mock_sce.rds" -o "$mock_data_dir/output_altExp.h5ad" -d
    if [ -f "$mock_data_dir/output_altExp.h5ad" ]; then
      echo "Test passed: altExp flag"
    else
      echo "Test failed: altExp flag"
    fi
    ;;
  "No")
    echo "Running no input test"
    Rscript -e "convert2anndata::cli_convert()" -o "$mock_data_dir/output_no_input.h5ad"
    if [ $? -ne 0 ]; then
      echo "Test passed: No input"
    else
      echo "Test failed: No input"
    fi
    ;;
  *)
    echo "Unknown input type"
    ;;
esac
