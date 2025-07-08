outpath <- "ArchR_TAM_anno"
setwd(outpath)

library(ArchR)
library(magrittr)
library(tidyverse)
addArchRGenome("hg38")
addArchRThreads(threads = 60)
library(BSgenome.Hsapiens.UCSC.hg38)

projHeme3 <- readRDS("/mnt/disk1/yu_new/project_3/Save-ProjHeme3/Save-ArchR-Project.rds")
idx <- which(projHeme3$anno == "Marcophages")
proj <- projHeme3[idx, ]


proj <- addIterativeLSI(proj,
    useMatrix = "TileMatrix", name = "IterativeLSI_TAM", iterations = 2,
    clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
    varFeatures = 25000, dimsToUse = 1:30
)

proj <- addClusters(proj, reducedDims = "IterativeLSI_TAM", method = "Seurat", resolution = 1.5)

proj <- addUMAP(proj, reducedDims = "IterativeLSI_TAM", name = "UMAP_TAM", nNeighbors = 30, minDist = 0.5)

seRNA <- readRDS("Macro_anno.rds")
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI_TAM",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "predicted.id",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


ggplot(data.frame(predictedScore_Un = proj$predictedScore_Un), aes(x = predictedScore_Un)) +
    geom_histogram(fill = "lightgrey", color = "black") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkred") +
    labs(title = "Efficiency for label transfer", x = "Predicted_Score", y = "Frequency") +
    theme_minimal()

proj <- addImputeWeights(proj)

)

plotEmbedding(proj,
    colorBy = "cellColData", name = "macro_anno2",
    embedding = "UMAP_TAM", size = 1, plotAs = "points", rastr = FALSE
)

proj <- addGroupCoverages(proj, sampleLabels = "Sample", groupBy = "macro_anno", force = TRUE)
proj <- addReproduciblePeakSet(proj, groupBy = "macro_anno", peakMethod = "Macs2", peaksPerCell = 500, force = TRUE)
projHeme_TAM3 <- addPeakMatrix(proj, force = TRUE)

markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme_TAM3,
    useMatrix = "PeakMatrix",
    groupBy = "macro_anno",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

peak_data <- as_tibble(getPeakSet(projHeme_TAM3))
peak_data$celltype <- sapply(strsplit(peak_data$GroupReplicate, "._."), "[", 1)


proj <- addMotifAnnotations(proj, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(markerTest, proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
motifsDo <- peakAnnoEnrichment(markerTest, proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")

df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[, 1]) %>%
    arrange(desc(mlog10Padj)) %>%
    mutate(rank = row_number())

pdf("motif_up.pdf", width = 5, height = 5)
ggplot(df_up, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(data = df_up[1:10, ], aes(label = TF), size = 1.5, color = "black") +
    theme_ArchR() +
    labs(x = "Rank Sorted TFs Enriched", y = "-log10(P-adj)")
dev.off()

df_down <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[, 1]) %>%
    arrange(desc(mlog10Padj)) %>%
    mutate(rank = row_number())

pdf("motif_down_SPP1.pdf", width = 5, height = 5)
ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(data = df_down[1:10, ], aes(label = TF), size = 1.5, color = "black") +
    theme_ArchR() +
    labs(x = "Rank Sorted TFs Enriched", y = "-log10(P-adj)")
dev.off()

proj <- addPeak2GeneLinks(proj, reducedDims = "IterativeLSI", useMatrix = "GeneIntegrationMatrix")
