library(SCENIC)
library(Seurat)
library(dplyr)
library(Matrix)
library(GENIE3)
# /storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/cardec_seurat/SCENIC/SCENIC_R
setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
data <- readRDS("RNA_avg_size20.rds")
motifs <- readRDS("signac.motifs.rds")
expr <- data@assays$RNA@data
glist <- rownames(expr)
clean_glist <- function(glist){
  idx <- grep("^RPS|XIST|RP11|^MT|^LINC|\\.|\\-", glist)
  glist <- glist[-idx]
  return(glist)
}
glist <- clean_glist(glist)
expr <- as.matrix(expr[glist, ])
glist <- rownames(expr)
portion <- rowSums(expr == 0)/ncol(expr)
glist <- glist[portion < 0.60]

exprMat_filtered <- expr[glist, ]

# get motif names
motif.names <- toupper(unlist(motifs@motif.names))
common.tf <- intersect(motif.names, rownames(exprMat_filtered))

setwd("/storage/chenlab/Users/qingnanl/signacatac/out/runGenie3/")

weightMat <- GENIE3(exprMat_filtered, regulators=common.tf, nCores = 20, verbose = T)

saveRDS(weightMat, "weightMat.rds")



