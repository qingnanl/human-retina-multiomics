#conda activate archr
#R
library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("hg19")
#inputFiles <- getTutorialData("Hematopoiesis")#check format
inputFiles<-c(D19D013_Fovea="~/19D013_fovea/outs/fragments.tsv.gz",
              D19D013_Macular="~/19D013_Macular/outs/fragments.tsv.gz",
              D19D014_Fovea="~/19D014_fovea/outs/fragments.tsv.gz",
              .....)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 2000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP", 
    LSIMethod = 1
)

archr <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "archR_test",
  copyArrows = TRUE 
)

archr <- filterDoublets(archr, filterRatio = 1.2)
saveArchRProject(ArchRProj = archr, outputDirectory = "Save-archr", load = FALSE)

archr <- addIterativeLSI(
    ArchRProj = archr,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 4,
    clusterParams = list( 
        resolution = c(0.5),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 15000,
    dimsToUse = 1:30
)
archr <- addUMAP(
    ArchRProj = archr,
    reducedDims = "IterativeLSI",
    name = "UMAP-LSI",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)
archr <- addClusters(
    input = archr,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.5,
    force = TRUE
)

p1 <- plotEmbedding(ArchRProj = archr, colorBy = "cellColData", name = "Sample", size = 0.001, plotAs = "points", labelMeans = FALSE,embedding = "UMAP-LSI")
p2 <- plotEmbedding(ArchRProj = archr, colorBy = "cellColData", name = "Clusters", size = 0.001, plotAs = "points", embedding = "UMAP-LSI")
plotPDF(p1,p2, name = "Plot-UMAP-LSI-Sample-Clusters.pdf", ArchRProj = archr, addDOC = FALSE, width = 8, height = 8)

saveArchRProject(ArchRProj = archr, outputDirectory = "Save-archr", load = FALSE)

# manual annotation based on markers 


# add annotation and remove doublets
library(dplyr)
major_type<-archr$Clusters
major_type<-recode(major_type, "C1" = "Rod", "C2" = "Rod", "C3" = "Doublet", "C4" = "Cone", "C5" = "MG-1", "C6" = "MG-2", "C7" = "Astrocyte", "C8" = "NN",
                      "C9" = "BC", "C10" = "BC", "C11" = "BC", "C12" = "BC", "C13" = "BC", "C14" = "BC", "C15" = "Doublet", "C16" = "BC",
                     "C17" = "BC", "C18" = "AC", "C19" = "AC", "C20" = "AC", "C21" = "AC", "C22" = "AC", "C23" = "AC", "C24" = "RGC", "C25" = "HC")
archr$major_type<-major_type
p2 <- plotEmbedding(ArchRProj = archr, colorBy = "cellColData", name = "major_type", size = 0.001, plotAs = "points",labelMeans = FALSE, sampleCells = 50000, embedding = "UMAP-LSI")
plotPDF(p2, name = "Plot-UMAP-LSI-major_type-2.pdf", ArchRProj = archr, addDOC = FALSE, width = 8, height = 8)

#
saveArchRProject(ArchRProj = archr, outputDirectory = "before_filtering", load = FALSE)
idxPass <- which(archr$major_type != "Doublet")
cellsPass <- archr$cellNames[idxPass]
archr<-archr[cellsPass, ]
# 
idxPass <- which(archr$DoubletScore < quantile(archr$DoubletScore, 0.975))
cellsPass <- archr$cellNames[idxPass]
archr<-archr[cellsPass, ]
#call peaks
archr <- addGroupCoverages(ArchRProj = archr, groupBy = "Clusters", force = TRUE)
pathToMacs2 <-  "~/anaconda3/envs/snap_atac/bin/macs2"
archr <- addReproduciblePeakSet(
    ArchRProj = archr,
    groupBy = "Clusters",force = TRUE,
    method = "q",
    cutOff = 0.01,
    additionalParams = "--nomodel",
    pathToMacs2 = pathToMacs2
)
archr <- addPeakMatrix(archr, force = TRUE)
saveArchRProject(ArchRProj = archr, outputDirectory = "post_filtering_peak_calling", load = FALSE)

