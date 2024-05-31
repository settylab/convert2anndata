library(SingleCellExperiment)

# Create a mock SingleCellExperiment object
counts <- matrix(1:100, nrow = 10)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
saveRDS(sce, "tests/testdata/mock_sce.rds")

# Create a mock Seurat object if Seurat is available
if (requireNamespace("Seurat", quietly = TRUE)) {
  library(Seurat)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
  saveRDS(seurat_obj, "tests/testdata/mock_seurat.rds")
}
