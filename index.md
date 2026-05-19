Bidirectional conversion between AnnData (.h5ad) and either Seurat or
SingleCellExperiment objects. Handles split layers (Seurat), assays,
dimensional reductions (obsm \<-\> reducedDims), metadata (obs/var \<-\>
colData/rowData), layers, and alternative experiments, aiming for a
faithful roundtrip.
