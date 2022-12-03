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
clean_glist <- function(glist){
  idx <- grep("^RPS|XIST|RP11|^MT|^LINC|\\.|\\-", glist)
  glist <- glist[-idx]
  return(glist)
}

data <- readRDS("atac_full_avg.rds")
DefaultAssay(data) <- "RNA"
glist <- rownames(data)
portion <- rowSums(data@assays$RNA@data==0)/ncol(data)
glist <- glist[portion < 0.70]
glist <- clean_glist(glist)

DefaultAssay(data) <- "peak"
data <- RegionStats(data, genome = BSgenome.Hsapiens.UCSC.hg19)
data <- FindTopFeatures(data, min.cutoff = 30)
# Idents(data) <- "major_type"
# data <- subset(data, downsample = 8000)
data <- LinkPeaks(
      object = data,
      peak.assay = "peak",
      expression.assay = "RNA",
      genes.use = glist,
      pvalue_cutoff = 1,
      distance = 5e+05,
      min.distance = 1000,
      score_cutoff = 0,
      min.cells = 30
)
test <- Links(object = data)
#saveRDS(test, "peak_rna_transfer_link_ds.rds")
saveRDS(test, "peak_rna_transfer_full_avg_link.rds")

#saveRDS(data, "atac4_major_ds_8000.rds")
