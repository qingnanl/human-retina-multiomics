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


setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
#data <- readRDS("processed_atac2.rds")
#print(Assays(data))
#DefaultAssay(data) <- "RNA"
#data[['ATAC']] <- NULL
#data[['peak']] <- NULL
#saveRDS(data, "signac_gscore_sr.rds")
data <- readRDS("signac_gscore_sr.rds")

run_integration <- function(rna, atac, file_name){
  rna<-FindVariableFeatures(rna, selection.method = "dispersion",nfeatures = 5000)
  atac<-FindVariableFeatures(atac, selection.method = "dispersion",nfeatures = 5000)
  features.use<-intersect(VariableFeatures(rna), VariableFeatures(atac))
  print(length(features.use))
  rnatmp <- CreateSeuratObject(rna@assays$RNA@data[features.use, ], meta.data = rna@meta.data)
  atactmp <- CreateSeuratObject(atac@assays$RNA@data[features.use, ], meta.data = atac@meta.data)
  rnatmp <- ScaleData(rnatmp, do.scale = F,
         do.center = F, features = rownames(rnatmp), 
         vars.to.regress = "sample")
  rnatmp@assays$RNA@data <- rnatmp@assays$RNA@scale.data
  atactmp <- ScaleData(atactmp, do.scale = F,
         do.center = F, features = rownames(atactmp), 
         vars.to.regress = c("sample", "nCount_peak"))
  atactmp@assays$RNA@data <- atactmp@assays$RNA@scale.data

  transfer.anchors <- FindTransferAnchors(
      reference = rnatmp,
      query = atactmp,
      normalization.method = "LogNormalize",
      features = features.use,
      reference.assay = "RNA",
      query.assay = "RNA",
      reduction = "cca"
  )
  refdata_rna <- GetAssayData(rna, assay = "RNA", slot = "data")
  imputation_rna <- TransferData(anchorset = transfer.anchors, refdata = refdata_rna, 
                                 weight.reduction = "cca", dims = 1:20)
  #predicted.labels <- TransferData(
  #  anchorset = transfer.anchors,
  #  refdata = rnatmp$subtype,
  #  weight.reduction = "cca",
  #  dims = 1:20
  #)
  saveRDS(imputation_rna, paste0(file_name, "_transfer_rna.rds"))
}

Idents(data) <- "major_type"

HC_atac <- subset(x = data, idents = c("HC"))
HC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/HC_final.h5")
setwd("/storage/chenlab/Users/qingnanl/signacatac/integration/")
run_integration(rna = HC_rna, atac = HC_atac, file_name = "HC")

BC_atac <- subset(x = data, idents = c("BC"))
BC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/BC_final.h5")
run_integration(rna = BC_rna, atac = BC_atac, file_name = "BC")

AC_atac <- subset(x = data, idents = c("AC"))
AC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/AC_final.h5")
run_integration(rna = AC_rna, atac = AC_atac, file_name = "AC")

RGC_atac <- subset(x = data, idents = c("RGC"))
RGC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/RGC_final.h5")
run_integration(rna = RGC_rna, atac = RGC_atac, file_name = "RGC")

Cone_atac <- subset(x = data, idents = c("Cone"))
Cone_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/Cone_final.h5")
run_integration(rna = Cone_rna, atac = Cone_atac, file_name = "Cone")

Rod_atac <- subset(x = data, idents = c("Rod"))
Rod_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/Rod_final.h5")
Rod_rna <- subset(Rod_rna, downsample = 50000)
run_integration(rna = Rod_rna, atac = Rod_atac, file_name = "Rod")

NN_atac <- subset(x = data, idents = c("MG", "Astrocyte"))
NN_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/NN_final.h5")
run_integration(rna = NN_rna, atac = NN_atac, file_name = "NN")
