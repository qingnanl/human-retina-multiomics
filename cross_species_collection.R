library(dplyr)
library(Seurat)
library(dior)
library(Matrix.utils)
setwd("/storage/chenlab/Users/qingnanl/humanre/out/cspecies")

#chicken
H_chi_gene<-read.delim("Human_chicken.txt")
H_chi_gene<-H_chi_gene[!duplicated(H_chi_gene$Chicken.gene.name), ]
H_chi_gene<-H_chi_gene[!duplicated(H_chi_gene$Gene.name), ]
rownames(H_chi_gene)<-H_chi_gene$Chicken.gene.name
H_chi_gene[H_chi_gene==""]<-NA
H_chi_gene<-H_chi_gene[complete.cases(H_chi_gene),]
chi <- read.csv("Chick_retina_atlas_expression_matrix.csv", row.names = 1)
chi<-chi[!duplicated(row.names(chi)), ]
chi <- chi[-grep("xcc", row.names(chi)), ]
colnames(chi) <- gsub("[.]", "-", colnames(chi))
chi_common<-intersect(rownames(H_chi_gene), rownames(chi))
H_chi_gene<-H_chi_gene[chi_common, ]
chi<-chi[chi_common, ]
rownames(chi)<-as.character(H_chi_gene$Gene.name)

chi_sr <- CreateSeuratObject(chi)
chi_meta <- read.csv("Chick_retina_atlas_meta.csv")
rownames(chi_meta) <- chi_meta$NAME
chi_sr@meta.data <- chi_meta

Idents(chi_sr) <- "Cluster"
chi_ct <- unique(as.character(chi_sr@meta.data$Cluster))
CAC <- subset(chi_sr, idents = grep("AC-", chi_ct, value = T))
saveRDS(CAC, "CAC.rds")
CBC <- subset(chi_sr, idents = grep("BP-", chi_ct, value = T))
saveRDS(CBC, "CBC.rds")
CHC <- subset(chi_sr, idents = grep("HC-", chi_ct, value = T))
saveRDS(CHC, "CHC.rds")
CRGC <- subset(chi_sr, idents = grep("RGC-", chi_ct, value = T))
saveRDS(CRGC, "CRGC.rds")
saveRDS(chi_sr, "chicken.rds")

# library(Matrix.utils)
H_Mac_gene<-read.csv("Human_macaque.txt")
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name), ]
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name.1), ]
rownames(H_Mac_gene)<-H_Mac_gene$Gene.name.1
H_Mou_gene<-read.csv("Human_mouse.txt")
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name), ]
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name.1), ]
rownames(H_Mou_gene)<-H_Mou_gene$Gene.name.1

#BC
HBC <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/BC_final.h5")
HBC@meta.data$subtype <- recode(HBC@meta.data$subtype, "BC0" = "FMB", "BC1" = "IMB", "BC2" = "RB",
                                   "BC3" = "DB2", "BC" = "DB1", "BC5" = "DB5", "BC6" = "DB4a", "BC7" = "DB3b",
                                   "BC8" = "DB3a", "BC9" = "BB_GB1", "BC10" = "DB6", "BC11" = "DB4b", "BC12" = "OFFx", "BC13" = "BB_GB2")
dior::write_h5(data = HBC, file = paste0("/storage/chenlab/Users/qingnanl/humanre/out/subtypes/BC_final.h5"))
HBC_meta<-HBC@meta.data
HBC_meta$species<-"human"
HBC_meta$cell_identity<-paste0("human_", HBC_meta$subtype)
HBC_mat<-HBC@assays$RNA@data
rm(HBC)

MABC<-readRDS("MBC.rds")
MABC_meta<-MABC@meta.data
MABC_meta$species<-"monkey"
MABC_meta$cell_identity<-paste0("monkey_", MABC_meta$ident)
MABC_mat<-MABC@assays$RNA@data
#here goes the -p / -n issue
#grep('-p', rownames(MABC_mat), values = T)
rownames(MABC_mat)<-gsub("-[pn]", "", rownames(MABC_mat))
MABC_mat<-aggregate.Matrix(MABC_mat,row.names(MABC_mat))
#MABC_mat<-t(sapply(by(MABC_mat,rownames(MABC_mat),colSums),identity))
rm(MABC)
ma_common<-intersect(rownames(H_Mac_gene), rownames(MABC_mat))
H_Mac_gene<-H_Mac_gene[ma_common, ]
MABC_mat<-MABC_mat[ma_common, ]
rownames(MABC_mat)<-as.character(H_Mac_gene$Gene.name)


MOBC<-readRDS("mouse_BC.rds")
Idents(MOBC)<-"label"
MOBC<-subset(MOBC, ident = c("AC", "Cone", "MG", "Rod", "Doublets/Contaminants"), invert = T)
MOBC_meta<-MOBC@meta.data
MOBC_meta$species<-"mouse"
MOBC_meta$cell_identity<-paste0("mouse_", MOBC_meta$label)
MOBC_mat<-MOBC@assays$RNA@data
rm(MOBC)

mo_common<-intersect(rownames(H_Mou_gene), rownames(MOBC_mat))
H_Mou_gene<-H_Mou_gene[mo_common, ]
MOBC_mat<-MOBC_mat[mo_common, ]
rownames(MOBC_mat)<-as.character(H_Mou_gene$Gene.name)


CBC_meta<-CBC@meta.data
CBC_meta$species<-"chicken"
CBC_meta$cell_identity<-paste0("chicken_", CBC_meta$Cluster)
CBC_mat<-CBC@assays$RNA@data



BC_mat_common<-Reduce(intersect, list(rownames(HBC_mat),rownames(MABC_mat),rownames(MOBC_mat), rownames(CBC_mat)))
HBC_mat<-HBC_mat[BC_mat_common, ]
HBC_meta<-HBC_meta[, c("species", "cell_identity")]
HBC_sr<-CreateSeuratObject(HBC_mat, meta.data = HBC_meta)
MABC_mat<-MABC_mat[BC_mat_common, ]
MABC_meta<-MABC_meta[, c("species", "cell_identity")]
MABC_sr<-CreateSeuratObject(MABC_mat, meta.data = MABC_meta)
MOBC_mat<-MOBC_mat[BC_mat_common, ]
MOBC_meta<-MOBC_meta[, c("species", "cell_identity")]
MOBC_sr<-CreateSeuratObject(MOBC_mat, meta.data = MOBC_meta)

CBC_mat<-CBC_mat[BC_mat_common, ]
CBC_meta<-CBC_meta[, c("species", "cell_identity")]
CBC_sr<-CreateSeuratObject(CBC_mat, meta.data = CBC_meta)


BC_species<-merge(HBC_sr, y = c(MABC_sr, MOBC_sr, CBC_sr))
saveRDS(BC_species, "BC_species.rds")
rm(BC_species)
rm(HBC_mat)
rm(HBC_sr)
rm(MABC_mat)
rm(MABC_sr)
rm(MOBC_mat)
rm(MOBC_sr)
rm(CBC_mat)
rm(CBC_sr)

H_Mac_gene<-read.csv("Human_macaque.txt")
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name), ]
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name.1), ]
rownames(H_Mac_gene)<-H_Mac_gene$Gene.name.1
H_Mou_gene<-read.csv("Human_mouse.txt")
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name), ]
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name.1), ]
rownames(H_Mou_gene)<-H_Mou_gene$Gene.name.1

#AC
HAC <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/AC_final.h5")
HAC_meta<-HAC@meta.data
HAC_meta$species<-"human"
HAC_meta$cell_identity<-paste0("human_", HAC_meta$subtype)
HAC_mat<-HAC@assays$RNA@data
rm(HAC)

MAAC<-readRDS("MAC.rds")
MAAC_meta<-MAAC@meta.data
MAAC_meta$species<-"monkey"
MAAC_meta$cell_identity<-paste0("monkey_", MAAC_meta$ident)
MAAC_mat<-MAAC@assays$RNA@data
#here goes the -p / -n issue
#grep('-p', rownames(MAAC_mat), values = T)
rownames(MAAC_mat)<-gsub("-[pn]", "", rownames(MAAC_mat))
MAAC_mat<-aggregate.Matrix(MAAC_mat,row.names(MAAC_mat))
#MAAC_mat<-t(sapply(by(MAAC_mat,rownames(MAAC_mat),colSums),identity))
rm(MAAC)
ma_common<-intersect(rownames(H_Mac_gene), rownames(MAAC_mat))
H_Mac_gene<-H_Mac_gene[ma_common, ]
MAAC_mat<-MAAC_mat[ma_common, ]
rownames(MAAC_mat)<-as.character(H_Mac_gene$Gene.name)


MOAC<-readRDS("mouse_AC.rds")
MOAC_meta<-MOAC@meta.data
MOAC_meta$species<-"mouse"
MOAC_meta$cell_identity<-paste0("mouse_", MOAC_meta$Cluster)
MOAC_mat<-MOAC@assays$RNA@data
rm(MOAC)

mo_common<-intersect(rownames(H_Mou_gene), rownames(MOAC_mat))
H_Mou_gene<-H_Mou_gene[mo_common, ]
MOAC_mat<-MOAC_mat[mo_common, ]
rownames(MOAC_mat)<-as.character(H_Mou_gene$Gene.name)
CAC_meta<-CAC@meta.data
CAC_meta$species<-"chicken"
CAC_meta$cell_identity<-paste0("chicken_", CAC_meta$Cluster)
CAC_mat<-CAC@assays$RNA@data

AC_mat_common<-Reduce(intersect, list(rownames(HAC_mat),rownames(MAAC_mat),rownames(MOAC_mat), rownames(CAC_mat)))
HAC_mat<-HAC_mat[AC_mat_common, ]
HAC_meta<-HAC_meta[, c("species", "cell_identity")]
HAC_sr<-CreateSeuratObject(HAC_mat, meta.data = HAC_meta)
MAAC_mat<-MAAC_mat[AC_mat_common, ]
MAAC_meta<-MAAC_meta[, c("species", "cell_identity")]
MAAC_sr<-CreateSeuratObject(MAAC_mat, meta.data = MAAC_meta)
MOAC_mat<-MOAC_mat[AC_mat_common, ]
MOAC_meta<-MOAC_meta[, c("species", "cell_identity")]
MOAC_sr<-CreateSeuratObject(MOAC_mat, meta.data = MOAC_meta)
CAC_mat<-CAC_mat[AC_mat_common, ]
CAC_meta<-CAC_meta[, c("species", "cell_identity")]
CAC_sr<-CreateSeuratObject(CAC_mat, meta.data = CAC_meta)
AC_species<-merge(HAC_sr, y = c(MAAC_sr, MOAC_sr, CAC_sr))
saveRDS(AC_species, "AC_species.rds")
rm(AC_species)
rm(HAC_mat)
rm(HAC_sr)
rm(MAAC_mat)
rm(MAAC_sr)
rm(MOAC_mat)
rm(MOAC_sr)
rm(CAC_mat)
rm(CAC_sr)
#
H_Mac_gene<-read.csv("Human_macaque.txt")
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name), ]
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name.1), ]
rownames(H_Mac_gene)<-H_Mac_gene$Gene.name.1
H_Mou_gene<-read.csv("Human_mouse.txt")
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name), ]
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name.1), ]
rownames(H_Mou_gene)<-H_Mou_gene$Gene.name.1

#HC
HHC <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/HC_final.h5")
HHC_meta<-HHC@meta.data
HHC_meta$species<-"human"
HHC_meta$cell_identity<-paste0("human_", HHC_meta$subtype)
HHC_mat<-HHC@assays$RNA@data
rm(HHC)

MAHC<-readRDS("MHC.rds")
MAHC_meta<-MAHC@meta.data
MAHC_meta$species<-"monkey"
MAHC_meta$cell_identity<-paste0("monkey_", MAHC_meta$ident)
MAHC_mat<-MAHC@assays$RNA@data
#here goes the -p / -n issue
#grep('-p', rownames(MAHC_mat), values = T)
rownames(MAHC_mat)<-gsub("-[pn]", "", rownames(MAHC_mat))
MAHC_mat<-aggregate.Matrix(MAHC_mat,row.names(MAHC_mat))
#MAHC_mat<-t(sapply(by(MAHC_mat,rownames(MAHC_mat),colSums),identity))
rm(MAHC)
ma_common<-intersect(rownames(H_Mac_gene), rownames(MAHC_mat))
H_Mac_gene<-H_Mac_gene[ma_common, ]
MAHC_mat<-MAHC_mat[ma_common, ]
rownames(MAHC_mat)<-as.character(H_Mac_gene$Gene.name)


MOHC<-readRDS("mouse_data.rds")
MOHC<-subset(MOHC, label == "HC")

MOHC_meta<-MOHC@meta.data
MOHC_meta$species<-"mouse"
MOHC_meta$cell_identity<-"mouse_HC"
MOHC_mat<-MOHC@assays$RNA@data
rm(MOHC)

CHC_meta<-CHC@meta.data
CHC_meta$species<-"chicken"
CHC_meta$cell_identity<-paste0("chicken_", CHC_meta$Cluster)
CHC_mat<-CHC@assays$RNA@data
#mo_common<-intersect(rownames(H_Mou_gene), rownames(MOHC_mat))
#H_Mou_gene<-H_Mou_gene[mo_common, ]
#MOHC_mat<-MOHC_mat[mo_common, ]
#rownames(MOHC_mat)<-as.character(H_Mou_gene$Gene.name)

HC_mat_common<-Reduce(intersect, list(rownames(HHC_mat),rownames(MAHC_mat),rownames(MOHC_mat),rownames(CHC_mat)))
HHC_mat<-HHC_mat[HC_mat_common, ]
HHC_meta<-HHC_meta[, c("species", "cell_identity")]
HHC_sr<-CreateSeuratObject(HHC_mat, meta.data = HHC_meta)
MAHC_mat<-MAHC_mat[HC_mat_common, ]
MAHC_meta<-MAHC_meta[, c("species", "cell_identity")]
MAHC_sr<-CreateSeuratObject(MAHC_mat, meta.data = MAHC_meta)
MOHC_mat<-MOHC_mat[HC_mat_common, ]
MOHC_meta<-MOHC_meta[, c("species", "cell_identity")]
MOHC_sr<-CreateSeuratObject(MOHC_mat, meta.data = MOHC_meta)
CHC_mat<-CHC_mat[HC_mat_common, ]
CHC_meta<-CHC_meta[, c("species", "cell_identity")]
CHC_sr<-CreateSeuratObject(CHC_mat, meta.data = CHC_meta)


HC_species<-merge(HHC_sr, y = c(MAHC_sr, MOHC_sr, CHC_sr))
saveRDS(HC_species, "HC_species.rds")
rm(HC_species)
rm(HHC_mat)
rm(HHC_sr)
rm(MAHC_mat)
rm(MAHC_sr)
rm(MOHC_mat)
rm(MOHC_sr)
rm(CHC_mat)
rm(CHC_sr)
## RGC
zib <- read.delim("expression_file2.tsv", row.names = 1)
H_zib_gene<-read.delim("Human_zebrafish.txt")
H_zib_gene<-H_zib_gene[!duplicated(H_zib_gene$Zebrafish.gene.name), ]
H_zib_gene<-H_zib_gene[!duplicated(H_zib_gene$Gene.name), ]
rownames(H_zib_gene)<-H_zib_gene$Zebrafish.gene.name
H_zib_gene[H_zib_gene==""]<-NA
H_zib_gene<-H_zib_gene[complete.cases(H_zib_gene),]

zib<-zib[!duplicated(row.names(zib)), ]

colnames(zib) <- gsub("[.]", "-", colnames(zib))
zib_common<-intersect(rownames(H_zib_gene), rownames(zib))
H_zib_gene<-H_zib_gene[zib_common, ]
zib<-zib[zib_common, ]
rownames(zib)<-as.character(H_zib_gene$Gene.name)

zib_sr <- CreateSeuratObject(zib)
zib_meta <- read.delim("metadata_file.tsv")
rownames(zib_meta) <- zib_meta$NAME
zib_sr@meta.data <- zib_meta
zib_sr <- subset(zib_sr, Cluster == "adult")


#
H_Mac_gene<-read.csv("Human_macaque.txt")
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name), ]
H_Mac_gene<-H_Mac_gene[!duplicated(H_Mac_gene$Gene.name.1), ]
rownames(H_Mac_gene)<-H_Mac_gene$Gene.name.1
H_Mou_gene<-read.csv("Human_mouse.txt")
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name), ]
H_Mou_gene<-H_Mou_gene[!duplicated(H_Mou_gene$Gene.name.1), ]
rownames(H_Mou_gene)<-H_Mou_gene$Gene.name.1

HRGC <- dior::read_h5(file = "/storage/chenlab/Users/qingnanl/humanre/out/subtypes/RGC_final.h5")
HRGC_meta<-HRGC@meta.data
HRGC_meta$species<-"human"
HRGC_meta$cell_identity<-paste0("human_", HRGC_meta$subtype)
HRGC_mat<-HRGC@assays$RNA@data
rm(HRGC)

MARGC<-readRDS("MRGC.rds")
MARGC_meta<-MARGC@meta.data
MARGC_meta$species<-"monkey"
MARGC_meta$cell_identity<-paste0("monkey_", MARGC_meta$ident)
MARGC_mat<-MARGC@assays$RNA@data
rm(MARGC)
ma_common<-intersect(rownames(H_Mac_gene), rownames(MARGC_mat))
H_Mac_gene<-H_Mac_gene[ma_common, ]
MARGC_mat<-MARGC_mat[ma_common, ]
rownames(MARGC_mat)<-as.character(H_Mac_gene$Gene.name)


MORGC<-readRDS("mouse_RGC.rds")
MORGC_meta<-MORGC@meta.data
MORGC_meta$species<-"mouse"
MORGC_meta$cell_identity<-paste0("mouse_", MORGC_meta$NAME)
MORGC_mat<-MORGC@assays$RNA@data
rm(MORGC)
mo_common<-intersect(rownames(H_Mou_gene), rownames(MORGC_mat))
H_Mou_gene<-H_Mou_gene[mo_common, ]
MORGC_mat<-MORGC_mat[mo_common, ]
rownames(MORGC_mat)<-as.character(H_Mou_gene$Gene.name)

CRGC_meta<-CRGC@meta.data
CRGC_meta$species<-"chicken"
CRGC_meta$cell_identity<-paste0("chicken_", CRGC_meta$Cluster)
CRGC_mat<-CRGC@assays$RNA@data

ZRGC_meta<-zib_sr@meta.data
ZRGC_meta$species<-"zebrafish"
ZRGC_meta$cell_identity<-paste0("zebrafish_RGC", ZRGC_meta$Subcluster)
ZRGC_mat<-zib_sr@assays$RNA@data


RGC_mat_common<-Reduce(intersect, list(rownames(HRGC_mat),rownames(MARGC_mat),rownames(MORGC_mat), rownames(CRGC_mat), rownames(ZRGC_mat)))
HRGC_mat<-HRGC_mat[RGC_mat_common, ]
HRGC_meta<-HRGC_meta[, c("species", "cell_identity")]
HRGC_sr<-CreateSeuratObject(HRGC_mat, meta.data = HRGC_meta)
MARGC_mat<-MARGC_mat[RGC_mat_common, ]
MARGC_meta<-MARGC_meta[, c("species", "cell_identity")]
MARGC_sr<-CreateSeuratObject(MARGC_mat, meta.data = MARGC_meta)
MORGC_mat<-MORGC_mat[RGC_mat_common, ]
MORGC_meta<-MORGC_meta[, c("species", "cell_identity")]
MORGC_sr<-CreateSeuratObject(MORGC_mat, meta.data = MORGC_meta)
CRGC_mat<-CRGC_mat[RGC_mat_common, ]
CRGC_meta<-CRGC_meta[, c("species", "cell_identity")]
CRGC_sr<-CreateSeuratObject(CRGC_mat, meta.data = CRGC_meta)
ZRGC_mat<-ZRGC_mat[RGC_mat_common, ]
ZRGC_meta<-ZRGC_meta[, c("species", "cell_identity")]
ZRGC_sr<-CreateSeuratObject(ZRGC_mat, meta.data = ZRGC_meta)


RGC_species<-merge(HRGC_sr, y = c(MARGC_sr, MORGC_sr, CRGC_sr, ZRGC_sr))
saveRDS(RGC_species, "RGC_species.rds")
rm(RGC_species)
rm(HRGC_mat)
rm(HRGC_sr)
rm(MARGC_mat)
rm(MARGC_sr)
rm(MORGC_mat)
rm(MORGC_sr)
