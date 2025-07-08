suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(future)
    library(Matrix)
    library(stringr)
})

nThreads <- 8
plan("multicore", workers = nThreads)

outpath <- "./3.scRNA"
setwd(outpath)

file_path <- "data"
scdata <- Read10X(file_path)
scobj <- CreateSeuratObject(
    counts = scdata,
    project = "scRNA",
    min.cells = 3,
    min.features = 200
)

{
    scobj <- NormalizeData(scobj)

    scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
    scobj <- ScaleData(scobj, features = rownames(scobj))

    scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj), reduction.name = "pca")
    scobj <- RunUMAP(scobj, reduction = "pca", dims = 1:30, reduction.name = "umap")
    scobj <- RunTSNE(scobj, reduction = "pca", dims = 1:30, reduction.name = "tsne")
    scobj <- FindNeighbors(scobj, dims = 1:30)
    scobj <- FindClusters(scobj, resolution = 0.5)
}

{
    scobj <- RunHarmony(scobj, reduction = "pca", group.by.vars = "Sample", reduction.save = "harmony")
    scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
}

{
    sample_col <- c("#d51f26", "#208045", "#2f2d66", "#6e4b9e", "#fbcb0a", "#c06cab", "#d8a767", "#8a9fd1", "#de6c3e", "#9983bd", "#3bbca8", "#90d5e4")
    names(sample_col) <- c(
        "H70", "H77", "H62", "H58a", "H63", "H65", "H23", "H30", "H38",
        "1HT1", "4HT1", "2HT1"
    )
    p1 <- DimPlot(scobj, group.by = "Sample", reduction = "umap", cols = sample_col)
    p2 <- DimPlot(scobj, group.by = "Sample", reduction = "umap_harmony", cols = sample_col)

    DimPlot(scobj, group.by = "Type", reduction = "umap_harmony", label = T, label.size = 6) + NoLegend()
    DimPlot(scobj, group.by = "seurat_clusters", reduction = "umap_harmony", label = T, label.size = 6) + NoLegend()

    DotPlot(scobj, features = c("HP", "KRT8", "KRT18")) # Hepatocytes
    DotPlot(scobj, features = c("ACTA2", "RGS5", "TAGLN")) # CAF
    DotPlot(scobj, features = c("GNG11", "VWF", "ENG")) # Endothelial

    DotPlot(scobj, features = c("CD79A", "CD19", "MS4A1", "MZB1", "IGHA1", "IGHG3")) # B
    DotPlot(scobj, features = c("CD3D", "CD3E", "CD8B")) # T
    DotPlot(scobj, features = c("NKG7", "KLRF1", "FGFBP2")) # NK

    DotPlot(scobj, features = c("CD1E", "CD1C", "FCER1A")) # Dendritic
    DotPlot(scobj, features = c("S100A8", "S100A9", "FCN1")) # Monocyte
    DotPlot(scobj, features = c("CD163", "CD68", "SLC40A1")) # Macrophage
}

scobj$anno <- scobj$seurat_clusters
scobj$anno <- recode(scobj$anno,
    "0" = "T cell",
    "1" = "T cell",
    "2" = "Hepatocytes",
    "3" = "Hepatocytes",
    "4" = "Hepatocytes",
    "5" = "T cell",
    "6" = "T cell",
    "7" = "Dendritic",
    "8" = "Macrophage",
    "9" = "T cell",
    "10" = "Endothelial cell",
    "11" = "B cell",
    "12" = "T cell",
    "13" = "Monocyte/Macrophage",
    "14" = "CAFs",
    "15" = "T cell",
    "16" = "Hepatocytes",
    "17" = "B cell",
    "18" = "T cell",
    "19" = "Hepatocytes",
    "20" = "Hepatocytes",
    "21" = "Hepatocytes",
    "22" = "B cell",
    "23" = "CAFs",
    "24" = "T cell",
    "25" = "Hepatocytes"
)


{
    Hepatocytes <- c("HP", "KRT8", "KRT18")
    names(Hepatocytes) <- NULL # names(Hepatocytes) = rep('Hepatocytes',length(Hepatocytes))
    Endothelial <- c("GNG11", "VWF", "ENG") # names(Endothelial) = rep('Endothelial',length(Endothelial))
    Fibroblast <- c("ACTA2", "RGS5", "TAGLN") # names(Fibroblast) = rep('Fibroblast',length(Fibroblast))

    T_cell <- c("CD3D", "CD3E", "CD8B") # names(T) = rep('T',length(T))
    NK <- c("NKG7", "KLRF1", "FGFBP2") # names(NK) = rep('NK',length(NK))
    B <- c("CD79A", "MZB1", "IGHG3") # names(B) = rep('B/Plasma',length(B))
    # Myeloid = c('S100A9', 'CD68', 'C1QC', 'C1QA'); #names(Myeloid) = rep('Myeloid',length(Myeloid))

    Dendritic <- c("CD1E", "CD1C", "FCER1A") # ,'DUSP4');
    Monocyte <- c("S100A8", "S100A9", "FCN1") # ,"EREG");
    Macrophage <- c("CD163", "CD68", "SLC40A1") # ,"FOLR2");
    # Mast = c('TPSB2','KIT','CPA3'); #names(Mast) = rep('Mast',length(Mast))
    marker_all <- c(T_cell, NK, Monocyte, Macrophage, Hepatocytes, Fibroblast, Endothelial, Dendritic, B)


    pdf("markerdotplot.pdf", width = 7, height = 5)
    DotPlot(scobj, features = marker_all, group.by = "anno", col.min = 0, col.max = 3, dot.scale = 4) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
        labs(x = NULL, y = NULL) +
        guides(size = guide_legend("Percentage")) +
        scale_color_gradientn(colours = c("#0072b5", "#fefedf", "#bc3c29"))
    dev.off()
}
