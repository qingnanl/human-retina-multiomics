library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(dplyr)
library(harmony)
library(dior)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SRAVG)
setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")

data <- readRDS("processed_atac4.rds")
# data <- readRDS("atac4_major_ds_8000.rds")
DefaultAssay(data) <- "RNA"

data@assays$ATAC <- NULL
#
start_time <- Sys.time()
data_avg <- sravg(object = data, dr_key = 'lsi', dr_dims = 2:30, group_size = 10,
                             group_within = 'subtype',peak_assay = "peak", peak_slot = "data",gex_slot = 'data',
                             extra_meta = c('nCount_RNA', 'nFeature_RNA', 'nCount_peak', 'nFeature_peak'))
end_time <- Sys.time()
end_time - start_time
#
saveRDS(data_avg, "atac_full_avg.rds")
