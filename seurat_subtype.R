#conda activate archr R18
library(Seurat)
library(dplyr)
library(Matrix)

data<-readRDS("~/cardec_sr.rds")
Idents(data)<-"major_type"
##BC
BC_RNA<-subset(data, idents = "BC")
BC_RNA<-FindVariableFeatures(BC_RNA, nfeatures = 5000)
saveRDS(VariableFeatures(BC_RNA), "./subtype/BC_vargenes.rds")
Idents(BC_RNA)<-"subtype"
BC_markers<-FindAllMarkers(BC_RNA, test.use = "bimod", logfc.threshold = 0.25,min.pct = 0.2,only.pos = T)
BC_markers<-subset(BC_markers, p_val_adj<0.05)
write.csv(BC_markers, "./subtype/BC_markers.csv")
avg<-sapply(split(colnames(BC_RNA@assays$RNA@data), BC_RNA$subtype),
                         function(cells) rowMeans(BC_RNA@assays$RNA@data[,cells]))
write.csv(avg, "./subtype/BC_avg.csv")
saveRDS(BC_RNA, "./subtype/BC.rds")