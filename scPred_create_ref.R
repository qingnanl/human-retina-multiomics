#R3
library(Seurat)
library(scPred)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
ref<-readRDS("/storage/singlecell/qingnanl/human_10x_process/visualize_genes/data_reduction.rds")
ref.subset<-subset(ref, downsample = 10000, idents = levels(ref@active.ident)[c(1:5, 7:10)])
  
ref.sce<-as.SingleCellExperiment(ref.subset)
ref_metadata <- as.data.frame(colData(ref.sce))
table(ref_metadata$cell_type)
set.seed(1234)
i <- createDataPartition(ref.sce$cell_type, p = 0.70, list = FALSE)
ref_counts <- logcounts(ref.sce)
ref_cpm  <- apply(ref_counts, 2, function(x) (x/sum(x))*1000000)
train_data <- ref_cpm[, i]
test_data <- ref_cpm[, -i]
saveRDS(ref.sce, "ref.sce.dspl.rds")
train_info <- ref_metadata[i, , drop = FALSE]
test_info <- ref_metadata[-i, , drop = FALSE]
scp <- eigenDecompose(train_data, n = 15)
scPred::metadata(scp) <- train_info
scp <- getFeatureSpace(scp, pVar = "cell_type")
pdf("plotEigen.pdf", width = 8, height = 6)
plotEigen(scp, group = "cell_type")
dev.off()
scp <- trainModel(scp, seed = 66)
res <- getTrainPred(scp)
pdf("plotTrained.pdf", width = 8, height = 6)
plotTrainProbs(scp)
dev.off()
scp <- scPredict(scp, newData = test_data, threshold = 0.9)
getPredictions(scp)
scp@predMeta <- test_info
crossTab(scp, true = "cell_type")
saveRDS(scp, "scp_new.rds")