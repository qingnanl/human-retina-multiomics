library(Seurat)
library(ArchR)
library(ChIPpeakAnno)
library(ChIPseeker)
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)
addArchRThreads(threads = 8)
addArchRGenome("hg19")
archr<-loadArchRProject(path = "~", force = FALSE, showLogo = TRUE)

#correlation between all the pairs
#cA, you need to compute the CoAccessibility before
cA2 <- getCoAccessibility(
    ArchRProj = archr,
    corCutOff = -1,
    resolution = 1000,
    returnLoops = FALSE
)
id_promoter<-which(archr@peakSet$peakType == "Promoter")
cA_pr_2<-cA2[cA2$queryHits %in% id_promoter, ]
cA_pr_2$query_peak<-paste0(as.character(archr@peakSet[cA_pr_2$queryHits]@seqnames), ":", archr@peakSet[cA_pr_2$queryHits]@ranges@start, "-", (archr@peakSet[cA_pr_2$queryHits]@ranges@start +500))
cA_pr_2$query_type<-archr@peakSet[cA_pr_2$queryHits]$peakType
cA_pr_2$query_ng<-archr@peakSet[cA_pr_2$queryHits]$nearestGene
cA_pr_2$subject_peak<-paste0(as.character(archr@peakSet[cA_pr_2$subjectHits]@seqnames), ":", archr@peakSet[cA_pr_2$subjectHits]@ranges@start, "-", (archr@peakSet[cA_pr_2$subjectHits]@ranges@start +500))
cA_pr_2$subject_type<-archr@peakSet[cA_pr_2$subjectHits]$peakType
cA_pr_2$subject_ng<-archr@peakSet[cA_pr_2$subjectHits]$nearestGene
cA_pr_2$query_ng_distance<-archr@peakSet[cA_pr_2$queryHits]$distToTSS
cA_pr_2$query_subject_distance<-abs((as.numeric(gsub(".*:", "", gsub("\\-.*","",  cA_pr_2$subject_peak))) + as.numeric(gsub(".*-", "", cA_pr_2$subject_peak)))/2 - 
         (as.numeric(gsub(".*:", "", gsub("\\-.*","",  cA_pr_2$query_peak))) + as.numeric(gsub(".*-", "", cA_pr_2$query_peak)))/2)
cA_pr_2<-cA_pr_2[cA_pr_2$query_ng_distance < 1000, ]
cA_pr_2$subject_ng_distance<-archr@peakSet[cA_pr_2$subjectHits]$distToTSS
cA_pr_2<-cA_pr_2[cA_pr_2$subject_ng_distance > 1000, ]#exclude promoter peaks from subject peaks
cA_pr_2$paste<-paste0(cA_pr_2$subject_peak, "_", cA_pr_2$query_ng)
#gene_peak_pair
gene_peak<-read.csv("~/cor_table_sub_type_pearson_fdr.csv", row.names = 1)#
gene_peak_pr_2<-gene_peak[gene_peak$distance >1000, ]#
gene_peak_pr_2$paste<-paste0(gene_peak_pr_2$peak, "_", gene_peak_pr_2$gene)
#find union
all_union<-unique(intersect(cA_pr_2$paste, gene_peak_pr_2$paste))#
cA_df<-cA_pr_2[cA_pr_2$paste %in% all_union, ]
cA_df<-cA_df[!duplicated(cA_df$paste), ]
cA_df<-cA_df[order(cA_df$paste),]
gp_df<-gene_peak_pr_2[gene_peak_pr_2$paste %in% all_union, ]
gp_df<-gp_df[order(gp_df$paste),]

cor_df<-data.frame(cA = cA_df$correlation, gp = gp_df$estimate)
pdf('./overlapping_peaks/all_union_cor.pdf', width =6, height = 4)
ggplot(cor_df, aes(cA, gp)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+geom_smooth(method = lm, linetype="dashed",
             color="darkred", size = 2)+theme_classic()+
    theme(
      plot.title = element_text(size=12, colour = 'black'),
    axis.text=element_text(size=15, colour = 'black'))+ geom_vline(xintercept = 0.5, linetype="dotted", 
                color = "blue", size=2)+ geom_hline(yintercept = 0.5, linetype="dotted", 
                color = "blue", size=2)#+
  #geom_point(shape = '.')
dev.off()

#Then filter the two list by cor >0.5