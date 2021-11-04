#conda activate archr R18
library(Seurat)
library(dplyr)
library(Matrix)
library(DropletUtils)
sample_list<-c("D028_13", "D027_13", "D026_13", "D021_13", "D019_13", "D018_13", "D017_13", "D013_13", "D009_13","D005_13",
               "D030_13", "19_D019", "19_D011", "19_D010","19_D009", "19_D008",  "19_D007", "19_D006","19_D005", "19_D003")
head<-"~/output_1/"
tail<-"_sr.rds"
setwd("/storage/chen/data_share_folder/human_20_snRNAseq/output_1/")

for (sample in sample_list){
  address<-paste0(head, sample, tail)
  data<-readRDS(address)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")#get mitochondria percentage
  data <- subset(data, subset = percent.mt < quantile(data@meta.data$percent.mt, 0.95))
  Idents(data)<-"scPred"
  RGC<-SubsetData(data, ident.use = "RGC")
  data_output_dir_1 <- "/storage/chen/data_share_folder/human_20_snRNAseq/output_1/RGC/"
  data_output_dir<-paste0(data_output_dir_1, sample, "_RGC")
  write10xCounts(x = RGC@assays$RNA@counts, version = "3",path = data_output_dir)
  }
