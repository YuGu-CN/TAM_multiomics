library(Seurat)
library(tidyverse)
library(magrittr)
library(harmony)

outpath = 'scRNA_TAM'
scobj = readRDS('scRNA.rds')
scobj_macro = subset(scobj,idents = "Marcophages")

{
  scobj_macro <- NormalizeData(scobj_macro)
  
  scobj_macro <- FindVariableFeatures(scobj_macro, selection.method = "vst", nfeatures = 2000)
  scobj_macro <- ScaleData(scobj_macro, features = rownames(scobj_macro))
  
  scobj_macro <- RunPCA(scobj_macro, features = VariableFeatures(object = scobj_macro),reduction.name = "pca")
  scobj_macro <- RunUMAP(scobj_macro,reduction = "pca", dims = 1:30, reduction.name = "umap")
  scobj_macro <- RunTSNE(scobj_macro,reduction = "pca", dims = 1:30, reduction.name = "tsne")
  scobj_macro <- FindNeighbors(scobj_macro, dims = 1:30)
  scobj_macro <- FindClusters(scobj_macro, resolution = 0.5)
}

{
  scobj_macro <- RunHarmony(scobj_macro,reduction = "pca",group.by.vars = "Sample",reduction.save = "harmony")
  scobj_macro <- RunUMAP(scobj_macro, reduction = "harmony", dims = 1:30,reduction.name = "umap_harmony")
  scobj_macro <- RunTSNE(scobj_macro,reduction = "harmony", dims = 1:30, reduction.name = "tsne_harmony")
}
sample_col = c('#d51f26','#208045','#2f2d66','#6e4b9e','#fbcb0a','#c06cab','#d8a767','#8a9fd1','#de6c3e','#9983bd','#3bbca8','#90d5e4')
names(sample_col) = c("H70", "H77", "H62", "H58a", "H63", "H65", "H23", "H30", "H38", 
                      "1HT1", "4HT1", "2HT1")
DimPlot(scobj_macro,reduction = 'umap',group.by = 'Sample',cols = sample_col)+DimPlot(scobj_macro,reduction = 'umap_harmony',group.by = 'Sample',cols = sample_col)


ref_marker <- read.csv("marker_Nature.csv")
    Mph_vec <- marker$Cell.cluster %>%
    unique() %>%
    str_subset(., "^Mph")

    Mph_df <- marker %>% filter(Cell.cluster %in% Mph_vec)
    Mph_df <- Mph_df %>%
    group_by(Cell.cluster) %>%
    slice_head(n = 30)
Mph_marker_list <- split(Mph_df$Gene, Mph_df$Cell.cluster)


plot_list <- list()
    for (i in 1:length(Mph_marker_list)) {
        marker_score_name <- names(Mph_marker_list[i])
        marker_score <- AddModuleScore(
        object = scobj,
        features = Mph_marker_list[i],
        name = marker_score_name

        plot_list[[i]] <- ggplot(my_data, aes_string(x = "UMAP_1", y = "UMAP_2", colour = paste0(marker_score_name, "1"))) +
        geom_point(size = 0.01) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red") +
        theme_classic()
        )
    }
    # find marker
    DefaultAssay(scobj) = 'RNA'
    subcluster_marker = FindAllMarkers(scobj,only.pos = T)
    subcluster_marker = subcluster_marker %>% group_by(cluster) %>% arrange(-avg_log2FC,-pct.1) %>% slice_head(n=50) # intersect marker with ref

    scobj$final_macro_anno <- scobj$seurat_clusters
    scobj$final_macro_anno <- recode(scobj$final_macro_anno,
      "0" = "SPP1+ Macro",
      "1" = "SLC40A1+ Macro",
      "2" = "TREM2+ Macro",
      "3" = "CXCL9+ Macro",
      "4" = "EREG+ MDSCs",
      "5" = "LYZ+ MDSCs",
      "6" = "CLEC10A+ Macro",
      "7" = "EREG+ MDSCs",
      "8" = "FCGR3A+ MDSCs",
      "9" = "HSP+ Macro",
      "10" = "Undefined",
      "11" = "STMN1+ Macro",
      "12" = "LTB+ MDSCs",
      "13" = "CLEC10A+ Macro",
      "14" = "SPP1+ Macro",
      "15" = "Undefined",
      "16" = "Kupffer cell",
      "17" = "Undefined",
      "18" = "CXCL9+ Macro",
      "19" = "Undefined",
      "20" = "CXCL9+ Macro"
    )
    scobj$final_macro_anno <- factor(scobj$final_macro_anno, levels = c("Kupffer cell", "CXCL9+ Macro", "SPP1+ Macro", "TREM2+ Macro", "CLEC10A+ Macro", "SLC40A1+ Macro", "STMN1+ Macro", "HSP+ Macro", "EREG+ MDSCs", "FCGR3A+ MDSCs", "LTB+ MDSCs", "LYZ+ MDSCs", "Undefined"))

    macro_color_vec <- c("#f24013", "#b2ac4a", "#3ca55c", "#4facb7", "#1adaf7", "#057ddd", "#16a38d", "#60cf98", "#ee8e50", "#a16e5f", "#dcc75f", "#f4a318", "grey")
      names(macro_color_vec) <- c("Kupffer cell", "CXCL9+ Macro", "SPP1+ Macro", "TREM2+ Macro", "CLEC10A+ Macro", "SLC40A1+ Macro", "STMN1+ Macro", "HSP+ Macro", "EREG+ MDSCs", "FCGR3A+ MDSCs", "LTB+ MDSCs", "LYZ+ MDSCs", "Undefined")

    (DimPlot(scobj, group.by = "final_macro_anno", cols = macro_color_vec, label = T, raster = F, raster.dpi = c(200, 200)) + theme_sc(x.label = "UMAP 1", y.label = "UMAP 2") + NoLegend() + theme(title = element_blank()))

    DimPlot(scobj, group.by = "final_macro_anno", split.by = 'new_section',cols = macro_color_vec, label = T, raster = F, raster.dpi = c(200, 200)) + theme_sc(x.label = "UMAP 1", y.label = "UMAP 2") + NoLegend() + theme(title = element_blank())
