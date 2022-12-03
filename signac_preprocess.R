library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM


padir <- "/storage/chenlab/Users/qingnanl/clean_h5ad/twenty_atac_upload/"
dir_list <- c("19_D003", "19_D006", "19_D008", "19_D010", "19D013_fovea_19", 
              "19D014_19", "19_D019", "D009_13", "D017_13", "D019_13", "D026_13", 
              "D028_13", "19_D005", "19_D007", "19_D009", "19_D011", "19D013_macular_19",
              "19D016_19", "D005_13", "D013_13", "D018_13", "D021_13", "D027_13")
              
# get a common peak set
getpeak <- function(dir){
  directory <- paste0(padir, dir, "/peaks.bed")
  # print(directory)
  a <- read.table(file = directory, col.names = c("chr", "start", "end"))
  a <- makeGRangesFromDataFrame(a)
  return(a)
}

peaklst <- c()
for (i in 1:length(dir_list)){
  peaklst <- append(peaklst, getpeak(dir_list[i]))
}

combined.peaks <- reduce(x = peaklst)
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 8000 & peakwidths > 100]

setwd("/storage/chenlab/Users/qingnanl/signacatac/out")
saveRDS(combined.peaks, "initial_combined_peaks.rds")

# count and create object
countpeak <- function(dir){
  directory1 <- paste0(padir, dir, "/singlecell.csv")
  directory2 <- paste0(padir, dir, "/fragments.tsv.gz")
  # print(directory)
  a <- read.table(file = directory1, 
                  stringsAsFactors = FALSE,
                  sep = ",",
                  header = TRUE,
                  row.names = 1)[-1, ]
  a <- a[a$passed_filters > 1000, ]
  b <- CreateFragmentObject(path = directory2,
                            cells = rownames(a))
  c <- FeatureMatrix(
  fragments = b,
  features = combined.peaks,
  cells = rownames(a))
  
  d <- CreateChromatinAssay(c, fragments = b)
  e <- CreateSeuratObject(d, assay = "ATAC", meta.data=a)
  e$dataset <- dir
  return(e)
}
obj_list <- list()
for (i in 1:length(dir_list)){
  obj_list[i] <- countpeak(dir_list[i])
  print(i)
}

combined <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(dir_list)],
  add.cell.ids = dir_list
)
saveRDS(combined, "initial_combined_atac.rds")
