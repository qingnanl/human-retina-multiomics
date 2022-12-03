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
  refdata_rna <- GetAssayData(rnatmp, assay = "RNA", slot = "data")
  imputation_rna <- TransferData(anchorset = transfer.anchors, refdata = refdata_rna, 
                                 weight.reduction = "cca", dims = 1:20)
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = rnatmp$subtype,
    weight.reduction = "cca",
    dims = 1:20
  )
  sub_rna<-CreateSeuratObject(refdata_rna)
  sub_rna@meta.data$tech<-"RNA"
  sub_rna@meta.data$sample <- rna$sample
  sub_rna@meta.data$cell_type<-rna$subtype
  sub_rna@meta.data$predicted.id <- rna$subtype
  sub_atac<-CreateSeuratObject(imputation_rna@data)
  sub_atac@meta.data$tech<-"ATAC"
  sub_atac@meta.data$cell_type<-"ATAC"
  sub_atac@meta.data$sample <- atac$sample
  sub_atac <- AddMetaData(object = sub_atac, metadata = predicted.labels)
  
  coembed <- merge(x = sub_rna, y = sub_atac)
  coembed <- ScaleData(coembed, features = rownames(coembed))
  coembed <- RunPCA(coembed, features = rownames(coembed))
  coembed<-FindNeighbors(coembed, reduction = "pca", dims = 1:20)
  coembed<-FindClusters(coembed, resolution = 0.5)
  coembed <- RunUMAP(coembed, dims = 1:10, reduction = 'pca')
  pdf(paste0(file_name, "_tech.pdf"), height = 4, width = 10)
  p <- DimPlot(coembed, split.by = "tech", group.by = "tech", raster = T)
  print(p)
  dev.off()
  
  pdf(paste0(file_name, "_clusters.pdf"), width = 6, height = 5)
  p <- DimPlot(coembed, group.by = "seurat_clusters", label = T, raster = T)
  print(p)
  dev.off()
  
  pdf(paste0(file_name, "_cell_type.pdf"), width = 6, height = 5)
  p <- DimPlot(coembed, group.by = "cell_type", label = T, repel = T, raster = T)
  print(p)
  dev.off()
  
  pdf(paste0(file_name, "_cell_type_nolab.pdf"), width = 6, height = 5)
  p <- DimPlot(coembed, group.by = "cell_type", label = F, raster = T)
  print(p)
  dev.off()
  
  pdf(paste0(file_name, "_tech_cell_type.pdf"), height = 4, width = 10)
  p <- DimPlot(coembed, split.by = "tech", group.by = "predicted.id", label = T, repel = T, raster = T) + NoLegend()
  print(p)
  dev.off()
  
  pdf(paste0(file_name, "_tech_cell_type_nolab.pdf"), height = 4, width = 15)
  p <- DimPlot(coembed, split.by = "tech", group.by = "predicted.id", label = F, repel = F, raster = T)
  print(p)
  dev.off()
  
  rna.count<-table(coembed$tech)[2]
  confusion.test.1<-table(cluster = coembed@meta.data$seurat_clusters[1:rna.count], label = coembed@meta.data$cell_type[1:rna.count])
  confusion.test.1<-data.frame(confusion.test.1)
  confusion.test.1<-group_by(confusion.test.1, cluster) %>% mutate(percent = Freq/sum(Freq))
  pdf(paste0(file_name, "_cluster_rna_plot.pdf"), width = 10)
  p <- ggplot(data =  confusion.test.1, aes(x = cluster, y = label, color = percent, size =percent)) +
    geom_point()+scale_color_gradient(low="blue", high="red")+scale_size_area()+
    labs( x = 'cluster',
          y = 'label')+theme_bw()
  print(p)
  dev.off()
  
  saveRDS(coembed, paste0(file_name, "_coemb.rds"))
}

Idents(data) <- "major_type"
#HC_atac <- subset(x = data, idents = c("HC"))
#HC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/HC_final.h5")
setwd("/storage/chenlab/Users/qingnanl/signacatac/integration/")
#run_integration(rna = HC_rna, atac = HC_atac, file_name = "HC")

BC_atac <- subset(x = data, idents = c("BC"))
BC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/BC_final.h5")
run_integration(rna = BC_rna, atac = BC_atac, file_name = "BC")

AC_atac <- subset(x = data, idents = c("AC"))
AC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/AC_final.h5")
run_integration(rna = AC_rna, atac = AC_atac, file_name = "AC")

RGC_atac <- subset(x = data, idents = c("RGC"))
RGC_rna <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/RGC_final.h5")
run_integration(rna = RGC_rna, atac = RGC_atac, file_name = "RGC")



