library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(dplyr)

setwd("/storage/chenlab/Users/qingnanl/signacatac/out/")
data <- readRDS("processed_atac1.rds")
#
peaks <- readRDS("peaks.rds")

filterpeak <- function(set){
  a <- set[set$score > 100 & set$neg_log10qvalue_summit >5, ]
  return(a)
}
peaklst <- c()
for (i in 1:length(peaks)){
  peaklst <- append(peaklst, filterpeak(peaks[[i]]))
}

combined.peaks <- reduce(x = peaklst)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 2500 & peakwidths > 200]
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
saveRDS(combined.peaks, "call_peaks_combined.rds")

DefaultAssay(data) <- "ATAC"
macs2_counts <- FeatureMatrix(
  fragments = Fragments(data),
  features = combined.peaks,
  cells = colnames(data)
)

saveRDS(macs2_counts, "macs2_counts.rds")
# create a new assay using the MACS2 peak set and add it to the Seurat object
data[["peak"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(data),
  annotation = Annotation(data)
)
saveRDS(data, "processed_atac2.rds")
