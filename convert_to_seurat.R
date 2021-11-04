#archr env; R18

# convert scanpy output (h5ad) to seurat object
library(Seurat)
library(dplyr)
library(Matrix)
library(SeuratDisk)
Convert("~/data_annotated.h5ad", dest = "h5seurat", overwrite = FALSE)
data <- LoadH5Seurat("~/data_annotated.h5seurat", assays = "data")
saveRDS(data, "./cardec_sr.rds")