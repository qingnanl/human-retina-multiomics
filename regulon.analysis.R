library(SCENIC)
library(GENIE3)
library(Seurat)
library(dplyr)
library(Matrix)
library(reshape2)
library(cluster)
library(ape)
library(GSEABase)
library(AUCell)
library(pheatmap)
library(ggplot2)
library(cowplot)


# /storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/cardec_seurat/SCENIC/SCENIC_R
setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
data <- readRDS("RNA_avg_size20.rds")
data@meta.data$major_type <- recode(data@meta.data$subtype, "ON_MGC" = "RGC", "OFF_MGC" = "RGC", "OFF_PGC" = "RGC", "RGC8" = "RGC",
                                                            "RGC7" = "RGC", "ON_PGC" = "RGC", "ipRGC1" = "RGC", "ipRGC2" = "RGC", 
                                                            "RGC13" = "RGC", "RGC10" = "RGC",  "RGC11" = "RGC", "RGC9" = "RGC", "RGC12" = "RGC",
                                                            "AC24" = "AC", "AC2" = "AC", "AC8" = "AC",  "AC32" = "AC", "AC28" = "AC", 
                                                            "AC3" = "AC", "AC12" = "AC", "AC5" = "AC", "AC13" = "AC", "AC19" = "AC", 
                                                            "AC15" = "AC", "AC17" = "AC", "AC26" = "AC", "AC0" = "AC", "AC6" = "AC","AC30" = "AC",       
                                                            "AC1"  = "AC",  "AC27" = "AC", "AC7" = "AC", "AC33" = "AC", "AC16" = "AC", "AC29" = "AC", 
                                                            "AC4"  = "AC", "AC20"  = "AC",  "AC18" = "AC",  "AC9" = "AC", "AC10" = "AC", 
                                                            "AC14" = "AC",  "AC22"  = "AC",  "AC34" = "AC", "AC11" = "AC", "AC25"  = "AC",
                                                              "AC21" = "AC", "AC31"  = "AC", "AC35" = "AC", "AC23"  = "AC", 
                                                              "DB4a"  = "BC", "DB3a" = "BC", "FMB" = "BC", "RB" = "BC",  "DB1"  = "BC", "DB5" = "BC", 
                                                              "DB2" = "BC", "DB3b" = "BC", "BB_GB1" = "BC", "BB_GB2" = "BC", 
                                                              "IMB" = "BC", "DB6"  = "BC",  "DB4b" = "BC", 
                                                              "OFFx" = "BC",  "ML_Cone" = "Cone", "S_Cone" = "Cone", "HC0" = "HC",  "HC1" = "HC")
data <- RunUMAP(data, reduction = "scvi", dims = 1:30, n.neighbors = 10)
# load regulon
# mp <- readRDS("motifPositions.rds")
# rg1 <- readRDS("TF_DF_full_filter.rds")
rg4 <- readRDS("regulon.updated.rds")
rg4$TF <- as.character(rg4$TF)
rg4$Target <- as.character(rg4$Target)


# cor.mat <- cor(as.matrix(data@assays$RNA@data))
expr <- as.matrix(data@assays$RNA@data)
genes <- unique(union(rg4$TF, rg4$Target))
expr <- expr[genes, ]
cor.mat <- cor(t(expr))

rg4$cor <- 0
for (i in 1:nrow(rg4)){
  a <- rg4[i, "TF"]
  b <- rg4[i, "Target"]
  c <- cor.mat[a, b]
  rg4[i, "cor"] <- c
  if (i %% 1000 == 0){
    print(i)
    }
}
rg4 <- rg4[rg4$cor > 0, ]
saveRDS(rg4, "regulon.updated.2.rds")

cells_rankings <- AUCell_buildRankings(data@assays$RNA@data)

tf.tar.lst <- list()
for (tf in unique(rg4$TF)){
  tf.tar.lst[[tf]] <- rg4[rg4$TF == tf, ]$Target
}

x <- c()
for (tf in names(tf.tar.lst)){
  g <- GeneSet(tf.tar.lst[[tf]], setName = tf)
  # x[[tf]] <- g
  x <- append(x, g)
}
x <- GeneSetCollection(x)

cells_AUC <- AUCell_calcAUC(x, cells_rankings)
cells_AUC_df <- as.data.frame(cells_AUC@assays@data@listData$AUC)
cells_AUC_df <- cells_AUC_df[rowSums(cells_AUC_df) > 0, ]
# cells_AUC_df <- apply(cells_AUC_df, 1, function(x) {(x - min(x) / max(x) - min(x))})

setwd("/storage/chenlab/Users/qingnanl/signacatac/figures/")

target.per.tf <- as.data.frame(table(rg4$TF))
target.per.tf$cat <- "TF"
median(target.per.tf$Freq)

tf.per.target <- as.data.frame(table(rg4$Target))
tf.per.target$cat <- "Target"
median(tf.per.target$Freq)

tf.target.box <- rbind(target.per.tf, tf.per.target)

pdf("tf.metrics.1.pdf", width = 5, height = 2.5)
ggplot(target.per.tf, aes(x=Freq))+
  geom_histogram(color="darkblue", fill="lightblue") + 
  geom_vline(aes(xintercept=median(Freq)),
            color="black", linetype="dashed", size=1) + 
  theme_classic()
dev.off()

pdf("tf.metrics.2.pdf", width = 5, height = 2.5)
ggplot(tf.per.target, aes(x=Freq))+
  geom_histogram(color="darkblue", fill="lightblue") + 
  geom_vline(aes(xintercept=median(Freq)),
            color="black", linetype="dashed", size=1) + 
  theme_classic()
dev.off()


tf.retain <- as.character(target.per.tf[target.per.tf$Freq > 20, ]$Var1)

cells_AUC_df <- cells_AUC_df[tf.retain, ]
cells_AUC_df.p <- apply(cells_AUC_df, 1, function(x) {(x - min(x)) / (max(x) - min(x))})
cells_AUC_df.p <- as.data.frame(t(cells_AUC_df.p))
annotation <- data.frame(cell_type = data$major_type)
rownames(annotation) <- colnames(cells_AUC_df)
pdf("tf.activity.cell.annotation.pdf", width = 6, height = 10)
pheatmap(cells_AUC_df.p, annotation = annotation, 
         border_color = NA, #scale = "row",
         show_rownames = F, show_colnames = F)
dev.off()
#
tf.cor <- cor(t(cells_AUC_df))


cl <- hclust(as.dist(tf.cor))#

pdf("tf.activity.cor.pdf", width = 6, height = 6)
p <- pheatmap(as.data.frame(tf.cor), annotation = ctr, 
         border_color = NA, # scale = "row",
         #cluster_rows = cl, cluster_cols = cl,
         #treeheight_col = 0.9,
         cutree_cols = 10, cutree_rows = 10, 
         show_rownames = F, show_colnames = F)
p
dev.off()
ctr <- cutree(p$tree_row, k = 10)
ctr <- as.data.frame(ctr)
ctr$module <- paste0("module_", ctr$ctr)
ctr$ctr <- NULL
ctr$gene <- rownames(ctr)

lst <- list()
for (i in unique(ctr$module)){
  # print(i)
  # n <- paste0("module_", i)
  lst[[i]] <- ctr[ctr$module == i, "gene"]
}

identical(colnames(cells_AUC_df), colnames(data))
auc.avg <- as.data.frame(sapply(split(colnames(cells_AUC_df), data$major_type),
                         function(cells) rowMeans(cells_AUC_df[,cells])))
auc.avg$TF <- rownames(auc.avg)
auc.avg$module <- ctr$module
auc.df <- as.data.frame(melt(auc.avg))
auc.df$module <- factor(auc.df$module, levels = c(paste0("module_", 1:10)))
write.csv(auc.df, "auc.df.csv")

pdf("tf.modules.scores.pdf", width = 6, height = 6)
ggplot(auc.df, aes(x=variable, y=value, fill=module)) + 
    geom_boxplot(outlier.size=0.5) +
    coord_flip() + 
    # facet_wrap(~module, scale="free") + 
    theme_classic()
dev.off()

pdf("tf.modules.scores.pdf", width = 5, height = 6)
ggplot(auc.df, aes(factor(variable), value, fill = module)) +
        # geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        geom_boxplot(outlier.size = 0.5, fatten = 0.5, lwd = 0.5) + 
        scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
        facet_grid(rows = vars(module), scales = "free", switch = "y") +
        theme_cowplot(font_size = 12) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
        # ggtitle("Identity on x-axis") + 
        xlab("Identity") + ylab("Expression Level")
dev.off()


setwd("/storage/chenlab/Users/qingnanl/signacatac/figures/")
rg4 <- readRDS("regulon.updated.2.rds")
library(igraph)
rg5 <- rg4[rg4$Target %in% rg4$TF, ]
# rg5 <- rg5[rg5$weight > 0.01, ]
r <- rg5[, 1:2]
g <- graph_from_data_frame(r, directed = T)
g.outd <- degree(g, mode = c("out"))
write.csv(as.data.frame(g.outd), "g.outd.csv")

pdf("out.degree.pdf", width = 5, height = 3)
hist(g.outd, breaks = 30)
dev.off()

pdf("out.degree.network.pdf", width = 8, height = 8)
plot(g, 
     vertex.label = NA,
     edge.color = 'black',
     vertex.size = g.outd * 0.5 + 0.001,
     edge.arrow.size = 0.05,
     layout = layout.fruchterman.reingold)#layout_nicely(g))
dev.off()

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}
jac.tf <- list()
for (tf in names(V(g))[g.outd > 15]){#c("PAX6", "MEIS2", "MEF2C", "ZEB1", "PBX1", "OTX2", "NFIA", "HLF", "BACH2", "VSX2", "FOXN3", "TCF4", "RORA")){
  tf.tar <- rg5[rg5$TF == tf, "Target"]
  for (i in 1:10){
    module.genes <- lst[[paste0("module_", i)]]
    jac.tf[[tf]][i] <- jaccard(tf.tar, module.genes)
    }
}
jac.tf <- do.call(rbind, jac.tf)
colnames(jac.tf) <- paste0("module_", 1:10)

pdf("tf.module.jac.pdf", width = 6, height = 6)
pheatmap(jac.tf, #annotation = ctr, 
         border_color = NA, scale = "row",
         #cluster_rows = cl, cluster_cols = cl,
         #treeheight_col = 0.9,
         # cutree_cols = 10, cutree_rows = 10, 
         show_rownames = T, show_colnames = T)

dev.off()