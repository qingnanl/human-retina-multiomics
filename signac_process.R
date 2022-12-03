library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)

# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# seqlevelsStyle(annotations) <- 'UCSC'
setwd("/storage/chenlab/Users/qingnanl/signacatac/out")
annotations <- readRDS("annotations.rds")

setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
data <- readRDS("initial_combined_atac.rds")

Annotation(data) <- annotations

# compute nucleosome signal score per cell
data <- NucleosomeSignal(object = data)
# saveRDS(data, "initial_combined_atac.rds")
# compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

setwd("/storage/chenlab/Users/qingnanl/signacatac/figures/")
data$high.tss <- ifelse(data$TSS.enrichment > 2, 'High', 'Low')
pdf("tss_enrich_filter.pdf")
TSSPlot(data, group.by = 'high.tss') + NoLegend()
dev.off()

data$nucleosome_group <- ifelse(data$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("fragmentHist.pdf")
FragmentHistogram(object = data, group.by = 'nucleosome_group')
dev.off()

data <- subset(
  x = data,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
data
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q50')
data <- RunSVD(data)
setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
saveRDS(data, "processed_atac.rds")
