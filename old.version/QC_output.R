#conda archr; R18; Rscript10
# run QC on data; output as 10x format
args<-commandArgs(TRUE)
library(Seurat)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(DropletUtils)
data_dir_1<-"/storage/singlecell/qingnanl/human_sn_10x_storage_only/"
data_dir_2<-args[1]#e.g. 10x_Fovea_Nu_19D014/19D014_fovea_outs/19D014_fovea_filtered_feature_bc_matrix/
data_dir<-paste0(data_dir_1, data_dir_2)
print(data_dir)
data <- Read10X(data.dir = data_dir)
data <- CreateSeuratObject(counts = data, min.cells = 5)#loose initial conditions to retain most genes
print(dim(data))
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")#get mitochondria percentage
data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < quantile(data@meta.data$percent.mt, 0.95))
print(dim(data))
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:15)
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.5)
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK<-as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations_1 <- data@meta.data$seurat_clusters
homotypic.prop_1 <- modelHomotypic(annotations_1)           
nExp_poi_1 <- round(0.1*nrow(data@meta.data))  ## Assuming 10% doublet formation rate


## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi_1, reuse.pANN = FALSE, sct = FALSE)
dblt_label<-paste0("DF.classifications_0.25_", pK, "_", nExp_poi_1)
data@meta.data$dblt<-data@meta.data[[dblt_label]]
data <- subset(data, subset = dblt == "Singlet")
print(dim(data))
#output

data_output_dir_1 <- "/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/"
data_output_dir_2<-gsub("_filtered_feature_bc_matrix/", "", data_dir_2)
data_output_dir_2<-gsub(".*/","",data_output_dir_2)
data_output_dir<-paste0(data_output_dir_1, data_output_dir_2)
write10xCounts(x = data@assays$RNA@counts, version = "3",path = data_output_dir)
