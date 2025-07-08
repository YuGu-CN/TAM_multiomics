suppressPackageStartupMessages({
    library(ArchR)
    library(dplyr)
    library(tidyr)
    library(mclust)
})

source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/sample_metadata.R"))
source(paste0(scriptPath, "/cluster_labels.R"))

outpath <- "./3.scATAC"
setwd(outpath)

library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 60)
inputFiles <- list.files("./input_arrow", full.names = T)

ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    filterTSS = 3, # Dont set this too high because you can always increase later
    filterFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    cleanTmp = TRUE,
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10,
    knnMethod = "UMAP",
)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = TRUE
)

# major celltype
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list( # See Seurat::FindClusters
        resolution = 1,
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30
)

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)

proj <- addImputeWeights(proj)

markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.57")

markerGenes <- c(
    "HP", "KRT8", "KRT18", # Hepatocytes
    "GNG11", "VWF", "ENG", # Endothelial
    "ACTA2", "RGS5", "TAGLN", # Fibroblasts
    "CD3D", "CD3E", "CD8B", # T
    "NKG7", "KLRF1", "FGFBP2", # NK
    "CD79A", "MZB1", "IGHG3", # B
    "CD1E", "CD1C", "FCER1A", # Dendritic
    "S100A8", "S100A9", "FCN1", # Monocyte
    "CD163", "CD68", "SLC40A1", "C1QA", "C1QC", "C1QB", # Macrophage
)

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.00",
    labelMarkers = markerGenes,
    binaryClusterRows = TRUE,
    clusterCols = TRUE,
    transpose = FALSE
)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(hm, name = "GeneScores-Marker-Heatmap", width = 6, height = 10, ArchRProj = proj, addDOC = FALSE)

clustNames <- list(
    "C1" = "Endothelial",
    "C2" = "Fibroblasts",
    "C3" = "Hepatocytes",
    "C4" = "T",
    "C5" = "NK",
    "C6" = "Hepatocytes",
    "C7" = "B",
    "C8" = "Dendritic",
    "C9" = "Macrophage",
    "C10" = "Macrophage",
)
proj$NamedClust <- clustNames[proj$Clusters] %>% unlist()
