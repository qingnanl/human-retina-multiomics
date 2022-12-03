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
peaks <- CallPeaks(
  object = data,
  assay = "ATAC",
  group.by = "major_type",
  outdir = "/storage/chenlab/Users/qingnanl/signacatac/peaktmp/",
  fragment.tempdir = "/storage/chenlab/Users/qingnanl/signacatac/peaktmp/",
  cleanup = FALSE,
  combine.peaks = FALSE,
  macs2.path = "/storage/chen/home/qingnanl/anaconda3/envs/snap_atac/bin/macs2"
)

saveRDS(peaks, "peaks.rds")
