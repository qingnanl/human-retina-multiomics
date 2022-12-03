library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)

setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
data <- readRDS("processed_atac.rds")


setwd("/storage/chenlab/Users/qingnanl/signacatac/figures/")

#

data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3)
# pdf("atac_clustering.pdf")
# DimPlot(object = data, label = TRUE) + NoLegend()
# dev.off()

gene.activities <- GeneActivity(data)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
data[['RNA']] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

DefaultAssay(data) <- 'RNA'
feature_markers <- c("PDE6A", "GNAT1", "ARR3", "PDE6C", "TFAP2B", "GAD1", "SLC6A9", "TRPM1", "NETO1", "ONECUT1",
                    "ONECUT2", "SEPT4", "RBPMS", "NEFM", "RLBP1", "SLN", "C1QA", "C1QB", "CD74",
                    "GFAP", "SLC1A3", "RPE65")
pdf("atac_feature_marker_vln.pdf", width = 15, height = 10)
VlnPlot(
  object = data,
  features = feature_markers,
  stack = TRUE, flip = TRUE
)
dev.off()
setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
saveRDS(data, "processed_atac1.rds")
