# Load ArchR (and associated libraries)
library(ArchR)
library(Seurat)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

# ---------------- Figure 1

# Get additional functions, etc.:
scriptPath <- "code"
source(paste0(scriptPath, "/plop2g_vecing_config.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/archr_helpers.R"))

atacNamedClustCmap <- readRDS(paste0(scriptPath, "/scATAC_NamedClust_cmap.rds")) %>% unlist()
rnaNamedClustCmap <- readRDS(paste0(scriptPath, "/scRNA_NamedClust_cmap.rds")) %>% unlist()
source(paste0(scriptPath, "/cluster_labels.R"))

atac.NamedClust <- as.list(c(
    "Hepatocytes", "Endothelial cell", "Fibroblast", "T cell",
    "NK cell", "B cell", "Dendritic", "Marcophages"
))
names(atac.NamedClust) <- c(
    "Hepatocytes", "Endothelial cell", "Fibroblast", "T cell",
    "NK cell", "B cell", "Dendritic", "Marcophages"
)
rna.NamedClust <- atac.NamedClust

BroadClust <- atac.NamedClust



LNatacOrder <- unlist(atac.NamedClust)[NatacOrder]
LNrnaOrder <- unlist(rna.NamedClust)[NrnaOrder]
namedClustAspect <- 1.6
fineClustAspect <- 1.6

wd <- "/mnt/disk1/yu_new/project_3-2/out/figure_1"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

plotDir <- wd

# Load Genome Annotations
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38

rna_proj <- readRDS("scRNA_anno.rds")
atac_proj <- loadArchRProject("Save-ProjHeme3", force = TRUE)

rna_proj$Sample <- str_c("HCC-T-", rna_proj$Sample)
rna_proj$Sample <- recode(rna_proj$Sample, "HCC-T-H58a" = "HCC-T-H58")

rna_proj$F1_anno <- recode(rna_proj$F1_anno, "Monocyte" = "Marcophages")
rna_proj$F1_anno <- factor(rna_proj$F1_anno, levels = c(
    "Hepatocytes", "Endothelial cell", "Fibroblast", "T cell",
    "NK cell", "B cell", "Dendritic", "Marcophages"
))
atac_proj$F1_anno <- atac_proj$predictedGroup_Un
atac_proj$F1_anno <- factor(atac_proj$F1_anno, levels = c(
    "Hepatocytes", "Endothelial cell", "Fibroblast", "T cell",
    "NK cell", "B cell", "Dendritic", "Marcophages"
))
featureSets <- list(
    "Hepatocytes" = c("HP", "KRT8", "KRT18"), # names(Hepatocytes) = rep('Hepatocytes',length(Hepatocytes))
    "Endothelial cell" = c("GNG11", "VWF", "ENG"), # names(Endothelial) = rep('Endothelial',length(Endothelial))
    "Fibroblast" = c("ACTA2", "RGS5", "TAGLN"), # names(Fibroblast) = rep('Fibroblast',length(Fibroblast))

    "T cell" = c("CD3D", "CD3E", "CD8B"), # names(T) = rep('T',length(T))
    "NK cell" = c("NKG7", "KLRF1", "FGFBP2"), # names(NK) = rep('NK',length(NK))
    "B cell" = c("CD79A", "MZB1", "IGHG3"), # names(B) = rep('B/Plasma',length(B))

    "Dendritic" = c("CD1E", "CD1C", "FCER1A"), # ,'DUSP4');
    "Marcophages" = c("CD163", "CD68", "SLC40A1") # ,"FOLR2");
)

NatacOrder <- c(
    "Hepatocytes", "Endothelial cell", "Fibroblast", "T cell",
    "NK cell", "B cell", "Dendritic", "Marcophages"
)
NrnaOrder <- NatacOrder


# color map
sample_cmap <- c(
    "#da8a8b", "#e4b565", "#b01e1d", "#d67121", "#502986", "#005437",
    "#b8a0cb", "#50932c", "#b2db7b", "#a7c1da", "#415fa0", "#797778"
)
names(sample_cmap) <- c(
    "HCC-T-1HT1", "HCC-T-2HT1", "HCC-T-4HT1", "HCC-T-H23", "HCC-T-H30", "HCC-T-H38",
    "HCC-T-H58", "HCC-T-H62", "HCC-T-H63", "HCC-T-H65", "HCC-T-H70", "HCC-T-H77"
)
broadClustCmap <- c("#700c11", "#d52226", "#a3c8dc", "#4d2772", "#f37a20", "#5478a3", "#379838", "#acd386")
names(broadClustCmap) <- c("T cell", "NK cell", "B cell", "Hepatocytes", "Fibroblast", "Endothelial cell", "Marcophages", "Dendritic")

# Dot plot of RNA cluster markers
count_mat <- GetAssayData(object = rna_proj, slot = "counts")
avgPctMat <- avgAndPctExpressed(count_mat, rna_proj$F1_anno, feature_normalize = TRUE, min_pct = 0)

# Subset to genes we care about:
subGenes <- featureSets %>% do.call("c", .)
avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes, ]
avgPctMat <- avgPctMat[avgPctMat$grp %in% NrnaOrder, ]

# Assign labels
avgPctMat$grp <- unlist(rna.NamedClust)[as.character(avgPctMat$grp)]

# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0

long_df <- avgPctMat
row_col <- "feature"
col_col <- "grp"
val_col <- "avgExpr"

# Determine cluster and gene order:
wide_df <- unmelt(avgPctMat, row_col = "feature", col_col = "grp", val_col = "avgExpr")
wide_df <- prep2g_vecyOrderMat(wide_df[, LNrnaOrder], clusterCols = FALSE)

grp_order <- colnames(wide_df$mat)
# gene_order <- rownames(wide_df$mat) %>% rev() # Reverse this if planning on using plot vertically
gene_order <- c(
    "HP", "KRT18", "KRT8", "ENG", "GNG11", "VWF", "ACTA2",
    "RGS5", "TAGLN", "CD3D", "CD3E", "CD8B", "FGFBP2", "KLRF1", "NKG7",
    "CD79A", "IGHG3", "MZB1", "CD1C", "CD1E", "FCER1A", "CD163",
    "CD68", "SLC40A1"
) %>% rev()

pdf(paste0(plotDir, "/RNA_NamedClust_markers_dot_plot_scalp.pdf"), width = 6, height = 10)
dotPlot(avgPctMat,
    xcol = "grp", ycol = "feature", color_col = "avgExpr", size_col = "pctExpr",
    xorder = grp_order, yorder = gene_order, cmap = cmaps_BOR$sunrise, aspectRatio = namedClustAspect
)
dev.off()

# Dot plot of GeneScoreMatrix cluster markers
{
    GSM_se <- getMatrixFromProject(atac_proj, useMatrix = "GeneScoreMatrix")
    GSM_mat <- assays(GSM_se)$GeneScoreMatrix
    rownames(GSM_mat) <- rowData(GSM_se)$name

    avgPctMat <- avgAndPctExpressed(GSM_mat[, getCellNames(atac_proj)], atac_proj$F1_anno, feature_normalize = TRUE, min_pct = 0)

    # Subset to genes we care about:
    subGenes <- featureSets %>% do.call("c", .)
    avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes, ]
    avgPctMat <- avgPctMat[avgPctMat$grp %in% NatacOrder, ]

    # Assign labels
    avgPctMat$grp <- unlist(atac.NamedClust)[as.character(avgPctMat$grp)]

    # Threshold min pct
    avgPctMat$pctExpr[avgPctMat$pctExpr < 10] <- 0

    # Determine cluster and gene order:
    wide_df <- unmelt(avgPctMat, row_col = "feature", col_col = "grp", val_col = "avgExpr")

    # Keep same clustering as from RNA:
    wide_df <- wide_df[, LNatacOrder]
    grp_order <- colnames(wide_df)

    pdf(paste0(plotDir, "/GSM_NamedClust_markers_dot_plot.pdf"), width = 6, height = 9)
    dotPlot(avgPctMat,
        xcol = "grp", ycol = "feature", color_col = "avgExpr", size_col = "pctExpr",
        xorder = grp_order, yorder = gene_order, cmap = cmaps_BOR$horizonExtra, aspectRatio = namedClustAspect
    )
    dev.off()
}


#############################################################################
# Stacked Bar Plots of Sample by Cluster
#############################################################################
{
    broadOrder <- names(broadClustCmap)

    rnaBroadOrder <- broadOrder[broadOrder %in% unique(rna_proj$F1_anno)]
    atacBroadOrder <- broadOrder[broadOrder %in% unique(atac_proj$F1_anno)]
    LrnaBroadOrder <- unlist(BroadClust)[rnaBroadOrder]
    LatacBroadOrder <- unlist(BroadClust)[atacBroadOrder]

    barWidth <- 0.9

    # RNA

    # BroadClust by Sample barplot
    clustBySamp <- fractionXbyY(unlist(BroadClust)[rna_proj$F1_anno], rna_proj$Sample, add_total = TRUE, xname = "Cluster", yname = "Sample")
    clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels = c(LrnaBroadOrder, "total"), ordered = TRUE)
    pdf(paste0(plotDir, "/rna_sampleByBroadClustBarPlot.pdf"), height = 5, width = 6)
    stackedBarPlot(clustBySamp, cmap = sample_cmap, namedColors = TRUE, barwidth = barWidth)
    dev.off()

    # ATAC

    # BroadClust by Sample barplot
    clustBySamp <- fractionXbyY(unlist(BroadClust)[atac_proj$F1_anno], atac_proj$Sample, add_total = TRUE, xname = "Cluster", yname = "Sample")
    clustBySamp$Cluster <- factor(clustBySamp$Cluster, levels = c(LatacBroadOrder, "total"), ordered = TRUE)
    pdf(paste0(plotDir, "/atac_sampleByBroadClustBarPlot.pdf"), height = 5, width = 6)
    stackedBarPlot(clustBySamp, cmap = sample_cmap, namedColors = TRUE, barwidth = barWidth)
    dev.off()
}

# --------------------- Figure 2
{

    projHeme_TAM2 <- addGroupCoverages(ArchRProj = projHeme_TAM1, sampleLabels = "Sample", groupBy = "macro_anno", force = T)

    # assign macs2 path
    pathToMacs2 <- findMacs2()

    # call peak
    projHeme_TAM2 <- addReproduciblePeakSet(
        ArchRProj = projHeme_TAM2,
        peakMethod = "Macs2",
        groupBy = "macro_anno",
        peaksPerCell = 500,
        pathToMacs2 = pathToMacs2, force = T
    )


    # getPeakSet(projHeme4)

    projHeme_TAM3 <- addPeakMatrix(projHeme_TAM2, force = TRUE)
    # getAvailableMatrices(projHeme_TAM3)

    # ArchR 识别cell-type marker peak
    markersPeaks <- getMarkerFeatures(
        ArchRProj = projHeme_TAM3,
        useMatrix = "PeakMatrix",
        groupBy = "macro_anno",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    saveArchRProject(ArchRProj = projHeme_TAM3, outputDirectory = "Save-projHeme_TAM3", load = FALSE)
}

{
    celltype <- c(
        "CXCL9+ Macro", "SLC40A1+ Macro", "Undefined", "Kupffer cell",
        "SPP1+ Macro", "TREM2+ Macro", "CLEC10A+ Macro", "STMN1+ Macro"
    )
    forlook <- as.tibble(getPeakSet(projHeme_TAM3))
    forlook$celltype <- sapply(strsplit(forlook$GroupReplicate, "._."), "[", 1)

    plots_list <- list()
    for (i in 1:length(celltype)) {
        tmp_celltype <- celltype[i]
        forlook_tmp <- forlook %>% filter(celltype == tmp_celltype)
        df_summary <- forlook_tmp %>%
            group_by(peakType) %>%
            summarise(count = n()) %>%
            mutate(percent = count / sum(count) * 100)
        p <- ggplot(df_summary, aes(x = "", y = count, fill = peakType)) +
            geom_bar(width = 1, stat = "identity", color = "white") +
            coord_polar(theta = "y") +
            geom_text(aes(label = paste0(round(percent, 1), "%", "  (", count, ")")),
                position = position_stack(vjust = 0.5), size = 4
            ) +
            scale_fill_brewer(palep2g_vece = "Set2") +
            theme_void() +
            labs(
                title = tmp_celltype,
                x = NULL, y = NULL, fill = "peakType"
            ) +
            theme( # legend.text = element_text(size = 12)
                plot.title = element_text(size = 16, hjust = 0.5),
            )
        plots_list[[i]] <- p
    }
    pdf("./macro_peakType.pdf", width = 12, height = 6)
    grid_plot <- do.call(grid.arrange, c(plots_list, ncol = 4))
    dev.off()

    plots_list <- list()
    for (i in 1:length(celltype)) {
        tmp_celltype <- celltype[i]
        forlook_tmp <- forlook %>% filter(celltype == tmp_celltype & peakType == "Distal")
        tmp_data <- data.frame(
            disp2g_vecoTSS = forlook_tmp$disp2g_vecoTSS
        )

        p <- ggplot(tmp_data, aes(x = disp2g_vecoTSS)) +
            geom_histogram(fill = "lightgrey", color = "black") +
            # geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkred")+
            theme_minimal() +
            theme(
                panel.border = element_rect(colour = "black", fill = NA, size = 1),
                axis.text.x = element_text(size = 10),
                axis.text.y = element_text(size = 10),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(size = 18, hjust = 0.5)
            ) +
            labs(
                title = tmp_celltype,
                x = "distance", y = "Frenquency"
            )
        plots_list[[i]] <- p
    }
    pdf("./macro__disp2g_vecoTSS.pdf", width = 12, height = 6)
    grid_plot <- do.call(grid.arrange, c(plots_list, ncol = 4))
    dev.off()

    forlook2 <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
    for (i in 1:length(celltype)) {
        tmp_celltype <- celltype[i]
        forlook2_tmp <- forlook2[[tmp_celltype]] %>% as.data.frame()
        print(paste0(tmp_celltype, " ", nrow(forlook2_tmp)))
    }
}

{
    projHeme_TAM3 <- readRDS("Save-ArchR-Project.rds")
    markersPeaks <- getMarkerFeatures(
        ArchRProj = projHeme_TAM3,
        useMatrix = "PeakMatrix",
        groupBy = "macro_anno",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    heatmapPeaks <- markerHeatmap(
        seMarker = markersPeaks,
        cutOff = "FDR <= 0.05 & Log2FC >= 1",
        transpose = TRUE
    )

    markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
    PeakSet_df <- getPeakSet(projHeme_TAM3) %>% as.tibble()

    downstream_genelist <- list()
    for (i in names(markerList@listData)) {
        print(i)
        tmp_marker_peak_df <- markerList@listData[[i]] %>% as.data.frame()
        tmp_df <- data.frame(
            order = NA,
            gene = NA,
            peaktype = NA
        )
        for (j in 1:nrow(tmp_marker_peak_df)) {
            j_start <- tmp_marker_peak_df[j, "start"]
            tmp_df[j, 1] <- j
            tmp_df[j, 2] <- PeakSet_df %>%
                filter(start == j_start & idx == tmp_marker_peak_df[j, 2]) %>%
                .$nearestGene
            tmp_df[j, 3] <- PeakSet_df %>%
                filter(start == j_start & idx == tmp_marker_peak_df[j, 2]) %>%
                .$peakType
        }
        downstream_genelist[[i]] <- tmp_df
    }

    downstream_CLEC10A <- downstream_genelist[["CLEC10A+ Macro"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_CLEC10A <- downstream_CLEC10A[!duplicated(downstream_CLEC10A$gene), ]
    downstream_CXCL9 <- downstream_genelist[["CXCL9+ Macro"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_CXCL9 <- downstream_CXCL9[!duplicated(downstream_CXCL9$gene), ]
    downstream_Kupffer <- downstream_genelist[["Kupffer cell"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_Kupffer <- downstream_Kupffer[!duplicated(downstream_Kupffer$gene), ]
    downstream_SLC40A1 <- downstream_genelist[["SLC40A1+ Macro"]] # %>% filter(peaktype != 'Exonic' & peaktype != 'Promoter' & peaktype != 'Intronic')
    downstream_SLC40A1 <- downstream_SLC40A1[!duplicated(downstream_SLC40A1$gene), ]
    downstream_SPP1 <- downstream_genelist[["SPP1+ Macro"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_SPP1 <- downstream_SPP1[!duplicated(downstream_SPP1$gene), ]
    downstream_STMN1 <- downstream_genelist[["STMN1+ Macro"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_STMN1 <- downstream_STMN1[!duplicated(downstream_STMN1$gene), ]
    downstream_TREM2 <- downstream_genelist[["TREM2+ Macro"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_TREM2 <- downstream_TREM2[!duplicated(downstream_TREM2$gene), ]
    downstream_Undefined <- downstream_genelist[["Undefined"]] %>% filter(peaktype != "Exonic" & peaktype != "Promoter" & peaktype != "Intronic")
    downstream_Undefined <- downstream_Undefined[!duplicated(downstream_Undefined$gene), ]

    write.csv(downstream_CLEC10A, "downstream_CLEC10A.csv", row.names = F)
    write.csv(downstream_CXCL9, "downstream_CXCL9.csv", row.names = F)
    write.csv(downstream_Kupffer, "downstream_Kupffer.csv", row.names = F)
    write.csv(downstream_SLC40A1, "downstream_SLC40A1.csv", row.names = F)
    write.csv(downstream_SPP1, "downstream_SPP1.csv", row.names = F)
    write.csv(downstream_STMN1, "downstream_STMN1.csv", row.names = F)
    write.csv(downstream_TREM2, "downstream_TREM2.csv", row.names = F)
    write.csv(downstream_Undefined, "downstream_Undefined.csv", row.names = F)

    intersect(downstream_SPP1$gene, downstream_CXCL9$gene) %>% as.data.frame()

    pdf("Specific peak heatmap.pdf")
    draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    dev.off()
}

veen_df <- data.frame(gene = allPeaksGR$nearestGene, peaktype = allPeaksGR$peakType)
veen_p <- veen_df %>%
    filter(peaktype == "Promoter") %>%
    dplyr::select(gene) %>%
    unlist() %>%
    unique() %>%
    as.data.frame() %>%
    na.omit()
veen_i <- veen_df %>%
    filter(peaktype == "Intronic") %>%
    dplyr::select(gene) %>%
    unlist() %>%
    unique() %>%
    as.data.frame() %>%
    na.omit()
veen_d <- veen_df %>%
    filter(peaktype == "Distal") %>%
    dplyr::select(gene) %>%
    unlist() %>%
    unique() %>%
    as.data.frame() %>%
    na.omit()
library(xlsx)
write.xlsx(veen_p, file = "intersect_3_element_gene.xlsx", sheetName = "veen_p", row.names = FALSE)
write.xlsx(veen_i, file = "intersect_3_element_gene.xlsx", append = TRUE, sheetName = "veen_i", row.names = FALSE)
write.xlsx(veen_d, file = "intersect_3_element_gene.xlsx", append = TRUE, sheetName = "veen_d", row.names = FALSE)

# -------------------------------------------- Figure 3
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
atac_proj <- loadArchRProject("/TAM_F3", force = TRUE)


corrCutoff <- 0.4
varCutoffATAC <- 0.25
varCutoffRNA <- 0.25
coAccCorrCutoff <- 0.4

allPeaksGR <- getPeakSet(atac_proj)
allPeaksGR$peakName <- (allPeaksGR %>%
    {
        paste0(seqnames(.), "_", start(.), "_", end(.))
    })
names(allPeaksGR) <- allPeaksGR$peakName

atac_proj <- addPeak2GeneLinks(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI"
)
atac_proj <- addCoAccessibility(
    ArchRProj = atac_proj,
    reducedDims = "IterativeLSI"
)

full_p2gGR <- as(peak2gene_list, "GRangesList") %>% unlist()
full_coaccessibility <- as(coaccessibility_list, "GRangesList") %>% unlist()

{
    p2g_vec <- data.frame(peak_name = p2gMat@listData[["Peak2GeneLinks"]]@listData$peak %>% unique())
    p2g_vec$peak_name <- str_replace_all(p2g_vec$peak_name, ":", "_")
    p2g_vec$peak_name <- str_replace_all(p2g_vec$peak_name, "-", "_")
    p2g_vec_ref <- as.data.frame(allPeaksGR)
    p2g_vec <- left_join(p2g_vec, p2g_vec_ref, by = c("peak_name" = "peakName"))
    # table(p2g_vec$peakType) %>% prop.table()
    p2g_vec_polar <- data.frame(
        Category = c("Distal", "Exonic", "Intronic", "Promoter"),
        Percentage = c(0.3355260, 0.1038689, 0.4059879, 0.1546172)
    )
    pdf("polar.pdf", height = 5, width = 5)
    ggplot(p2g_vec_polar, aes(x = 4, y = Percentage, fill = Category)) +
        geom_col() +
        coord_polar(theta = "y") +
        xlim(c(0.2, 4.5)) +
        theme_void() +
        scale_fill_brewer(palep2g_vece = "Pastel1") +
        labs(fill = "Category") +
        theme(legend.title = element_blank())
    dev.off()

    p2g_vec_bed <- p2g_vec$peak_name %>% as.data.frame()
    colnames(p2g_vec_bed) <- "peak_name"
    p2g_vec_bed <- p2g_vec_bed %>% separate(peak_name, into = c("chr", "start", "end"), sep = "_")
    write.table(p2g_vec_bed[, c(1, 2, 3)], "p2glink_peak.bed", sep = "\t", quote = F, col.names = F, row.names = F)
    rm(p2g_vec)
    rm(p2g_vec_ref)
    rm(p2g_vec_polar)
    rm(p2g_vec_bed)

    # linux: cd Veen_bed
    # linux:intervene venn -i ./*.bed --output ./
    # linux:bedtools intersect -a p2glink_peak.bed -b ENCODE.bed FANTOM5.bed SCenhancer.bed -wa | sort -u > res.bed
    # 49783/65342
}

{
    matches <- getMatches(atac_proj, "Motif")
    r1 <- SummarizedExperiment::rowRanges(matches)
    rownames(matches) <- paste(seqnames(r1), start(r1), end(r1), sep = "_")
    matches <- matches[names(allPeaksGR)]

    clusters <- unique(kclust_df$kclust) %>% sort()

    enrichList <- lapply(clusters, function(x) {
        cPeaks <- kclust_df[kclust_df$kclust == x, ]$peakName %>% unique()
        ArchR:::.computeEnrichment(matches, which(names(allPeaksGR) %in% cPeaks), seq_len(nrow(matches)))
    }) %>% SimpleList()
    names(enrichList) <- clusters
}

{
    library(BSgenome.Hsapiens.UCSC.hg38)
    atac_proj <- addMotifAnnotations(ArchRProj = atac_proj, motifSet = "cisbp", name = "Motif")
}

#--------------------------- Figure 4
library(Scissor)
library(Seurat)
library(Augur)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(survival)

{
    scobj_scissor <- scobj
    rm(list = setdiff(ls(), c("exp_lihc", "pheno_lihc", "scobj_scissor")))


    scobj_scissor <- subset(scobj_scissor, final_macro_anno %in% c(
        "SPP1+ Macro", "CXCL9+ Macro", "CLEC10A+ Macro", "HSP+ Macro",
        "SLC40A1+ Macro", "STMN1+ Macro", "TREM2+ Macro"
    ))
    scobj_scissor <- subset(scobj_scissor, new_section == "Normal", invert = T)
    Idents(scobj_scissor) <- "final_macro_anno"
    if (identical(colnames(exp_lihc), pheno_lihc$ID)) {
        bulk_dataset <- exp_lihc
        phenotype <- pheno_lihc[, c("OS.time", "OS")]
        colnames(phenotype) <- c("time", "status")
        sc_dataset <- scobj_scissor
    }
    DefaultAssay(sc_dataset) <- "RNA"
    sc_dataset <- FindNeighbors(sc_dataset, dims = 1:30)
    sc_dataset@graphs$RNA_snn <- sc_dataset@graphs$integrated_snn
    infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype,
        alpha = NULL,
        family = "cox", Save_file = "Scissor_HCC_survival.RData"
    )

    Scissor_select <- rep("Background cells", ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)
    Scissor_select[infos1$Scissor_pos] <- "Scissor+ cell"
    Scissor_select[infos1$Scissor_neg] <- "Scissor- cell"
    sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
    prop.table(table(sc_dataset$final_macro_anno, sc_dataset$scissor), margin = 1)
    sc_dataset$scissor <- factor(sc_dataset$scissor, levels = c("Scissor+ cell", "Scissor- cell", "Background cells"))
    DimPlot(sc_dataset,
        reduction = "umap", group.by = "scissor", cols = c("#b2182b", "#2166ac", "lightgrey"), pt.size = 0.6,
        shuffle = T
    )
}

{
    ## construct reference matrix
    load("scobj.rdata")
    scobj_ciber <- scobj

    exp_matrix <- scobj_ciber@assays$RNA@data
    metadata <- scobj_ciber@meta.data
    identical(row.names(scobj_ciber@meta.data), colnames(exp_matrix))

    scobj_ciber$new_major_anno <- factor(scobj_ciber$new_major_anno, levels = dput(unique(scobj_ciber$new_major_anno)))
    cellnames <- as.character(levels(scobj_ciber$new_major_anno))
    sample_name <- c()
    nsample <- 500

    for (i in 1:length(cellnames)) {
        set.seed(1)
        newmetadata <- metadata[metadata$new_major_anno == cellnames[i], ]
        sample_name <- c(
            sample_name,
            row.names(newmetadata[sample(1:nrow(newmetadata), nsample, replace = F), ])
        )
        sample_name <- gsub("\\.[0-9]", "", sample_name)
        rm(newmetadata)
    }

    cell_names <- metadata[sample_name, ]$new_major_anno
    exp_matrix <- exp_matrix[, sample_name]
    colnames(exp_matrix) <- cell_names
    exp_matrix <- as.matrix(exp_matrix)
    adjustdata <- function(data) {
        data <- cbind(rownames(data), data)
    }
    exp_matrix <- adjustdata(exp_matrix)
    colnames(exp_matrix)[1] <- "GeneSymbol"

    ## prepare bulk expression data,and used the online tool CIBERSORTx to evaluate the cell type proportions

    ## compare
    ciber_spp1_res <- read.csv(".CIBERSORTx_Job1_Results.csv")
    lihc_clinical <- read.csv(".TCGA cohort_clinical data.csv")

    ciber_spp1_res <- ciber_spp1_res[str_sub(ciber_spp1_res$Mixture, 14, 15) == "01", ]
    lihc_clinical <- lihc_clinical %>% dplyr::select(ID, Stage, time)
    ciber_spp1_res <- dplyr::left_join(ciber_spp1_res, lihc_clinical, by = c("Mixture" = "ID"))
    ciber_spp1_res$Stage <- dplyr::recode(ciber_spp1_res$Stage,
        "1" = "I",
        "2" = "II",
        "3" = "III/IV",
        "4" = "III/IV"
    )
    ciber_spp1_res <- ciber_spp1_res %>% dplyr::select(ID, SPP1..Macro, Stage, time, os.event)

    ggplot(ciber_spp1_res, aes(x = Stage, y = SPP1..Macro, fill = Stage)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.signif",
            comparisons = list(c("I", "II"), c("I", "III/IV"), c("II", "III/IV"))
        ) +
        labs(
            title = "Comparison of SPP1..Macro Across Different Stages",
            x = "Disease Stage",
            y = "SPP1..Macro"
        ) +
        theme_minimal() +
        scale_fill_manual(values = c("I" = "lightblue", "II" = "lightgreen", "III/IV" = "salmon"))

    ciber_spp1_res$SPP1_Group <- ifelse(
        ciber_spp1_res$SPP1..Macro > median(ciber_spp1_res$SPP1..Macro, na.rm = TRUE), "High", "Low"
    )
    surv_obj <- Surv(time = ciber_spp1_res$time, event = ciber_spp1_res$os.event)
    surv_fit <- survfit(surv_obj ~ SPP1_Group, data = ciber_spp1_res)
    surv_plot <- ggsurvplot(
        surv_fit,
        data = ciber_spp1_res,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        legend.title = "SPP1..Macro Group",
        legend.labs = c("Low", "High"),
        xlab = "Time (Months)",
        ylab = "Survival Probability",
        title = "Survival Analysis by SPP1..Macro Group"
    )
}

{
    df <- left_join(df_counts, df_total, by = "Region") %>%
        mutate(across(Kupffer_cell:SPP1_Macro, ~ .x / TotalCells * 100)) %>%
        select(-TotalCells)

    df_long <- df %>%
        pivot_longer(cols = -Region, names_to = "CellType", values_to = "Percent") %>%
        mutate(
            Region = factor(Region, levels = c("Normal", "TumorEdge", "TumorCore")),
            CellType = factor(CellType, levels = c("Kupffer_cell", "CXCL9_Macro", "SPP1_Macro"))
        )

    df_long$Region <- factor(df_long$Region, levels = c("Normal", "TumorEdge", "TumorCore"))
    df_long$CellType <- factor(df_long$CellType,
        levels = c("Kupffer_cell", "CXCL9_Macro", "SPP1_Macro")
    )

    ggline(df_long,
        x = "Region", y = "Percent", color = "CellType",
        add = c("mean_se"), # 加标准误
        size = 1.2, point.size = 3, linetype = "solid"
    ) +
        scale_color_manual(values = c(
            "Kupffer_cell" = "#1C91C0",
            "CXCL9_Macro" = "#999999",
            "SPP1_Macro" = "#D73027"
        )) +
        labs(y = "Proportion of cells (%)", x = NULL, title = "Cell Proportion by Region") +
        theme_minimal(base_size = 14) +
        theme(
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            legend.title = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA)
        )
}

{
    DefaultAssay(scobj) <- "RNA"
    marker_df <- FindMarkers(scobj, only.pos = F, ident.1 = "SPP1+ Macro", ident.2 = "CXCL9+ Macro", min.pct = 0.1, logfc.threshold = 0)

    {
        df_markers <- marker_df %>%
            rownames_to_column("gene") %>%
            mutate(
                pct_diff = 100 * (pct.1 - pct.2),
                direction = case_when(
                    avg_log2FC > 1.5 & pct_diff > 10 ~ "up_in_A",
                    avg_log2FC < -1 & pct_diff < -10 ~ "up_in_B",
                    TRUE ~ "nonsig"
                ),
                label = ifelse(direction != "nonsig", gene, NA)
            )
    }

    g <- ggplot(df_markers, aes(x = pct_diff, y = avg_log2FC)) +
        geom_point(aes(color = direction), alpha = 0.7, size = 1.5) +
        geom_label_repel(
            data = df_markers %>% filter(direction == "up_in_A"),
            aes(label = gene),
            fill = "#E6550D", color = "white", size = 3.5, fontface = "bold",
            box.padding = 0.4, segment.color = "grey60", max.overlaps = 30
        ) +
        geom_label_repel(
            data = df_markers %>% filter(direction == "up_in_B"),
            aes(label = gene),
            fill = "#3182BD", color = "white", size = 3.5, fontface = "bold",
            box.padding = 0.4, segment.color = "grey60", max.overlaps = 30
        ) +
        scale_color_manual(values = c("up_in_A" = "#E6550D", "up_in_B" = "#3182BD", "nonsig" = "grey80")) +
        labs(
            x = "Δ Percentage Difference",
            y = "Log-Fold Change",
            title = "Differential Expression Between Cell Types"
        ) +
        theme_classic(base_size = 14) +
        theme(
            legend.position = "none",
            axis.text = element_text(color = "black", size = 12),
            axis.title = element_text(size = 13)
        )
}
####################### Figure 5

{
    atac_proj <- addImputeWeights(atac_proj)

    trajectory_SPP1 <- c("Kupffer cell", "SPP1+ Macro")
    atac_proj <- addTrajectory(
        ArchRProj = atac_proj,
        name = "Tra_SPP1",
        groupBy = "macro_anno",
        trajectory = trajectory_SPP1,
        embedding = "UMAP_TAM",
        force = TRUE
    )
    plotTrajectory(atac_proj, trajectory = "Tra_SPP1", colorBy = "cellColData", name = "Tra_SPP1", embedding = "UMAP_TAM")[[1]]


    trajMM_SPP1 <- getTrajectory(atac_proj, name = "Tra_SPP1", useMatrix = "MotifMatrix", log2Norm = FALSE)
    trajGSM_SPP1 <- getTrajectory(atac_proj, name = "Tra_SPP1", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
    trajGIM_SPP1 <- getTrajectory(atac_proj, name = "Tra_SPP1", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)

    p1 <- plotTrajectoryHeatmap(trajMM_SPP1, pal = paletteContinuous("solarExtra"))
    p2 <- plotTrajectoryHeatmap(trajGSM_SPP1, pal = paletteContinuous("horizonExtra"))
    p3 <- plotTrajectoryHeatmap(trajGIM_SPP1, pal = paletteContinuous("blueYellow"))
    plotPDF(p1, p2, p3, name = "Plot-Kup-SPP1-Traj-Heatmaps.pdf", ArchRProj = atac_proj, width = 6, height = 8)


    corGIM_MM_SPP1 <- correlateTrajectories(trajGIM_SPP1, trajMM_SPP1, corCutOff = 0.4, varCutOff1 = 0.6, varCutOff2 = 0.6)
    idxToRemove <- grep("deviations", corGIM_MM_SPP1$correlatedMappings$name2)
    if (length(idxToRemove) > 0) {
        corGIM_MM_SPP1$correlatedMappings <- corGIM_MM_SPP1$correlatedMappings[-idxToRemove, ]
    }

    genes_SPP1 <- corGIM_MM_SPP1$correlatedMappings$name1
    motifs_SPP1 <- corGIM_MM_SPP1$correlatedMappings$name2
    trajGIM2_SPP1 <- trajGIM_SPP1[genes_SPP1, ]
    trajMM2_SPP1 <- trajMM_SPP1[motifs_SPP1, ]

    combinedMat_SPP1 <- plotTrajectoryHeatmap(trajGIM2_SPP1, returnMat = TRUE, varCutOff = 0)
    assay_combined_SPP1 <- t(scale(t(assay(trajGIM2_SPP1)))) + t(scale(t(assay(trajMM2_SPP1))))
    assay(trajGIM2_SPP1, withDimnames = FALSE) <- assay_combined_SPP1

    rowOrder_SPP1 <- match(rownames(combinedMat_SPP1), rownames(trajGIM2_SPP1))
    ht1 <- plotTrajectoryHeatmap(trajGIM2_SPP1, pal = paletteContinuous("blueYellow"), rowOrder = rowOrder_SPP1)
    ht2 <- plotTrajectoryHeatmap(trajMM2_SPP1, pal = paletteContinuous("solarExtra"), rowOrder = rowOrder_SPP1)

    plotTrajectory(atac_proj, trajectory = "Tra_SPP1", colorBy = "GeneScoreMatrix", name = "Phagocytic process", continuousSet = "horizonExtra")
}

{
    keep_clusters <- c("Kupffer cell", "CXCL9+ Macro", "SPP1+ Macro")
    seurat_dorothea <- subset(scobj, final_macro_anno %in% keep_clusters, invert = FALSE)
    seurat_dorothea$final_macro_anno <- factor(seurat_dorothea$final_macro_anno, levels = c("Kupffer cell", "CXCL9+ Macro", "SPP1+ Macro"))

    seurat_dorothea$type <- seurat_dorothea$final_macro_anno
    Idents(seurat_dorothea) <- "type"

    data("dorothea_hs", package = "dorothea")
    regulon <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

    seurat_dorothea <- run_viper(seurat_dorothea, regulon,
        options = list(
            method = "scale", minsize = 4,
            eset.filter = FALSE, verbose = FALSE
        )
    )

    DefaultAssay(seurat_dorothea) <- "dorothea"
    seurat_dorothea <- ScaleData(seurat_dorothea)

    viper_scores_df <- GetAssayData(seurat_dorothea, slot = "scale.data", assay = "dorothea") %>%
        as.data.frame(check.names = FALSE) %>%
        t()

    CellsClusters <- data.frame(
        cell = names(Idents(seurat_dorothea)),
        cell_type = as.character(Idents(seurat_dorothea)),
        check.names = FALSE
    )

    viper_scores_long <- viper_scores_df %>%
        as.data.frame() %>%
        rownames_to_column("cell") %>%
        gather(tf, activity, -cell) %>%
        inner_join(CellsClusters, by = "cell")

    summarized_viper_scores <- viper_scores_long %>%
        group_by(tf, cell_type) %>%
        summarise(avg = mean(activity), std = sd(activity), .groups = "drop")

    top10_per_celltype <- summarized_viper_scores %>%
        group_by(cell_type) %>%
        slice_max(order_by = avg, n = 20, with_ties = FALSE) %>%
        ungroup()
}

{
    kgrps <- unique(knn_groups$group)


    pmatPsB <- lapply(kgrps, function(x) {
        use_cells <- knn_groups[knn_groups$group == x, ]$cell_name
        Matrix::rowSums(pmat[, use_cells])
    }) %>% do.call(cbind, .)
    colnames(pmatPsB) <- kgrps


    gsmatPsB <- lapply(kgrps, function(x) {
        use_cells <- knn_groups[knn_groups$group == x, ]$cell_name
        Matrix::rowMeans(gsmat[, use_cells])
    }) %>% do.call(cbind, .)
    colnames(gsmatPsB) <- kgrps


    gimatPsB <- lapply(kgrps, function(x) {
        use_cells <- knn_groups[knn_groups$group == x, ]$cell_name
        Matrix::rowMeans(gimat[, use_cells])
    }) %>% do.call(cbind, .)
    colnames(gimatPsB) <- kgrps


    scaled_pmatPsB <- t(t(pmatPsB) / colSums(pmatPsB)) * 10000
    L2gsmatPsB <- sparseLogX(gsmatPsB, logtype = "log2", scale = FALSE) #
    L2gimatPsB <- sparseLogX(gimatPsB, logtype = "log2", scale = FALSE) #
    psb_labels <- sapply(kgrps, function(kg) {
        clustLabels <- ccd[knn_groups$cell_name[knn_groups$group == kg], ]$macro_anno
        names(sort(table(clustLabels), decreasing = TRUE))[1]
    })

    {
        zmat <- assay(seZ)


        motifDeltaMat <- matrix(NA, nrow = nrow(zmat), ncol = ncol(zmat))
        rownames(motifDeltaMat) <- rownames(zmat)
        colnames(motifDeltaMat) <- colnames(zmat)


        for (i in seq_len(ncol(zmat))) {
            ref <- zmat[, i]
            delta <- abs(zmat - ref)
            delta[, i] <- NA
            motifDeltaMat[, i] <- apply(delta, 1, max, na.rm = TRUE) #
        }
        motifDeltaMat <- as.data.frame(motifDeltaMat)
        motifDeltaMat$TF_name <- rowData(seZ)$name
    }


    {
        motifSE <- getMatrixFromProject(atac_proj, useMatrix = "MotifMatrix", binarize = FALSE)
        z_motif <- assays(motifSE)$z
        rownames(z_motif) <- rowData(motifSE)$name
        z_motif <- z_motif[, getCellNames(atac_proj)] # Need to get


        z_motifPsB <- lapply(kgrps, function(x) {
            use_cells <- knn_groups[knn_groups$group == x, ]$cell_name
            Matrix::rowMeans(z_motif[, use_cells]) # GImatrix is already scaled
        }) %>% do.call(cbind, .)
        colnames(z_motifPsB) <- kgrps
    }


    {
        L2gimatPsB <- L2gsmatPsB
    }


    {
        rownames(z_motifPsB) <- sub("_.*", "", rownames(z_motifPsB))
        common_rows <- intersect(rownames(z_motifPsB), rownames(L2gimatPsB))
        common_rows <- common_rows[order(match(common_rows, rownames(z_motifPsB)))]

        z_sub <- z_motifPsB[common_rows, , drop = FALSE]
        L2g_sub <- L2gimatPsB[common_rows, , drop = FALSE]

        spp1_cols <- names(psb_labels)[psb_labels == "SPP1+ Macro"]
        z_spp1 <- z_sub[, spp1_cols, drop = FALSE]
        # z_spp1 = log2(z_spp1+1)
        L2g_spp1 <- L2g_sub[, spp1_cols, drop = FALSE]
    }

    {
        rowCorTest_simple <- function(X, Y, padjMethod = "fdr", min = 5, use = "complete") {
            stopifnot(nrow(X) == nrow(Y))

            results <- lapply(seq_len(nrow(X)), function(i) {
                x <- X[i, ]
                y <- Y[i, ]

                if (sum(!is.na(x) & !is.na(y)) >= min) {
                    ct <- cor.test(x, y, use = use, method = "spearman")
                    return(c(cor = ct$estimate, pval = ct$p.value))
                } else {
                    return(c(cor = NA, pval = NA))
                }
            })

            df <- as.data.frame(do.call(rbind, results))
            df$padj <- p.adjust(df$pval, method = padjMethod)
            rownames(df) <- rownames(X)
            df$motif <- rownames(df)
            return(df)
        }
    }

    {
        rowCorTest_simple_df <- rowCorTest_simple(z_spp1, L2g_spp1)
    }

    {
        rownames(motifDeltaMat) <- sub("_.*", "", motifDeltaMat$TF_name)
        motifDeltaMat_sub <- motifDeltaMat[common_rows, , drop = FALSE]
        rowCorTest_simple_df$Zscore <- motifDeltaMat_sub$`SPP1+ Macro`
    }

    {
        library(ggplot2)
        library(dplyr)
        library(ggrepel)

        zscore_threshold <- quantile(rowCorTest_simple_df$Zscore, 0.75, na.rm = TRUE)


        rowCorTest_simple_df <- rowCorTest_simple_df %>%
            mutate(TFRegulator = ifelse(cor.rho > 0.25 & Zscore > zscore_threshold, "YES", "NO"))



        zscore_threshold <- quantile(rowCorTest_simple_df$Zscore, 0.75, na.rm = TRUE)
        rowCorTest_simple_df <- rowCorTest_simple_df %>%
            mutate(TFRegulator = ifelse(cor.rho > 0.25 & Zscore > zscore_threshold, "YES", "NO"))


        top_tf <- rowCorTest_simple_df %>% filter(TFRegulator == "YES")


        pdf("SPP1_PTF.pdf", width = 5, height = 5)
        ggplot(rowCorTest_simple_df, aes(x = cor.rho, y = Zscore)) +
            geom_point(aes(color = TFRegulator), size = ifelse(rowCorTest_simple_df$TFRegulator == "YES", 2.5, 1.5), alpha = 0.7) +
            scale_color_manual(values = c("NO" = "#A9A9A9", "YES" = "#D73027")) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.6) +
            ggrepel::geom_text_repel(
                data = top_tf,
                aes(label = motif),
                size = 3.5,
                box.padding = 0.3,
                max.overlaps = 30,
                min.segment.length = 0,
                color = "#D73027",
                segment.color = "grey60"
            ) +
            labs(
                x = "Correlation to gene expression",
                y = "TF Motif Delta",
                color = "Positive Regulator"
            ) +
            theme_classic(base_size = 14) +
            theme(
                legend.position = "none",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                axis.text = element_text(color = "black"),
                axis.title = element_text(size = 13)
            )
        dev.off()
    }


    {
        kupffer_cols <- names(psb_labels)[psb_labels == "Kupffer cell"]
        z_kupffer <- z_sub[, kupffer_cols, drop = FALSE]
        L2g_kupffer <- L2g_sub[, kupffer_cols, drop = FALSE]


        rowCorTest_kupffer_df <- rowCorTest_simple(z_kupffer, L2g_kupffer)


        rownames(motifDeltaMat) <- sub("_.*", "", motifDeltaMat$TF_name)
        motifDeltaMat_sub <- motifDeltaMat[common_rows, , drop = FALSE]
        rowCorTest_kupffer_df$Zscore <- motifDeltaMat_sub$`Kupffer cell`


        zscore_threshold_kupffer <- quantile(rowCorTest_kupffer_df$Zscore, 0.75, na.rm = TRUE)
        rowCorTest_kupffer_df <- rowCorTest_kupffer_df %>%
            mutate(TFRegulator = ifelse(cor.rho > 0.5 & Zscore > zscore_threshold_kupffer, "YES", "NO")) %>%
            filter(Zscore <= 5)

        top_tf_kupffer <- rowCorTest_kupffer_df %>% filter(TFRegulator == "YES")


        library(ggplot2)
        library(dplyr)
        library(ggrepel)

        pdf("Kupffer_PTF.pdf", width = 5, height = 5)
        ggplot(rowCorTest_kupffer_df, aes(x = cor.rho, y = Zscore)) +
            geom_point(aes(color = TFRegulator), size = ifelse(rowCorTest_kupffer_df$TFRegulator == "YES", 2.5, 1.5), alpha = 0.7) +
            scale_color_manual(values = c("NO" = "#A9A9A9", "YES" = "#D73027")) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.6) +
            ggrepel::geom_text_repel(
                data = top_tf_kupffer,
                aes(label = motif),
                size = 3.5,
                box.padding = 0.3,
                max.overlaps = 30,
                min.segment.length = 0,
                color = "#D73027",
                segment.color = "grey60"
            ) +
            labs(
                x = "Correlation to gene expression",
                y = "TF Motif Delta",
                color = "Positive Regulator"
            ) +
            theme_classic(base_size = 14) +
            theme(
                legend.position = "none",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 11),
                axis.text = element_text(color = "black"),
                axis.title = element_text(size = 13)
            )
        dev.off()
    }
}

# ------------------------------------ Figure 6
{
library(Seurat)
library(Matrix)
library(data.table)
library(tidyverse)

data_dir <- "cellranger_out"
folder_dir <- file.path(list.dirs(data_dir, recursive = F), "outs/raw_peak_bc_matrix")

macro_path = read.csv('macro_barcode.csv',header = T,row.names = 1)
macro_index = rownames(macro_path)

macro_list = list()
for (i in 1:length(folder_dir)) {
    # print(list.files(folder_dir[i]))
    sample_id = str_split(folder_dir[i],'/',simplify = T) %>% .[,8]
    mtx_path <- file.path(folder_dir[i], list.files(folder_dir[i])[2])
    cl_path <- file.path(folder_dir[i], list.files(folder_dir[i])[1])
    rl_path <- file.path(folder_dir[i], list.files(folder_dir[i])[3])
    mtx <- readMM(mtx_path)
    cl <- read.table(cl_path,
        header = F
    )[, 1]
    rl <- read.table(rl_path,
        header = F
    )
    rl$V4 <- stringr::str_c(rl$V1, rl$V2, rl$V3, sep = "_")
    rl <- rl$V4
    rownames(mtx) <- rl
    colnames(mtx) <- cl
    colnames(mtx) <- str_c(sample_id,'#',colnames(mtx))

    macro_mtx = mtx[,which(colnames(mtx) %in% macro_index)]
    macro_list[[sample_id]] = macro_mtx
}

macro_mtx =  do.call(cbind.data.frame,macro_list)
}

## ```python
{
    from APEC import clustering,plot,generate,convert
    EX_FILE = 'matrix_2mtx/APEC'
    convert.convert_10X('matrix_2mtx/macro_peak_bc_matrix', EX_FILE)
    clustering.build_accesson(EX_FILE, ngroup=600)
    clustering.cluster_byAccesson(EX_FILE, nc=0, norm='probability')
    plot.plot_tsne(EX_FILE, rs=0)
    plot.correlation(EX_FILE, cell_label='cluster', clip=[0,1])
    generate.gene_score(EX_FILE, genome_gtf='hg38_RefSeq_genes.gtf', distal=20000)
    generate.get_nearby_genes(EX_FILE)
    generate.differential_feature(EX_FILE, feature='accesson', target='0', vs='all')
    generate.super_enhancer(EX_FILE, super _range=500000, p_cutoff=0.05)
}

{
    SE_grange = read.csv('APEC/result/potential_super_enhancer.csv',sep = '\t')
    SE_grange_fil = SE_grange %>% filter(p.value<0.01)
SE_grange_fil_3 = SE_grange_fil %>% separate(position, into = c("chr", "start", "end"), sep = "[:\\-]")
SE_grange_fil_3$start <- as.numeric(SE_grange_fil_3$start)
SE_grange_fil_3$end <- as.numeric(SE_grange_fil_3$end)
SE_gr <- GRanges(
  seqnames = SE_grange_fil_3$chr,
  ranges = IRanges(start = SE_grange_fil_3$start, end = SE_grange_fil_3$end)
)

    SE_signal_df = data.frame(matrix(ncol = length(colnames(tmp_signal_matrix)), nrow = 0));colnames(SE_signal_df) = colnames(signal_groupSE_mtx)

for (i in 1:length(SE_gr)) {
  tmp_SE_gr = SE_gr[i]
  tmp_SE_name = paste0(c('SE',seqnames(tmp_SE_gr) %>% as.character(),ranges(tmp_SE_gr) %>% as.character()),collapse = '_')
  enhancer_hit = findOverlaps(fragment_info_gr,tmp_SE_gr) %>% queryHits(.)
  tmp_enhancer_name = fragment_info[enhancer_hit,'fragname']
  tmp_signal_matrix = signal_groupSE_mtx[tmp_enhancer_name,]
  tmp_column_means <- colMeans(tmp_signal_matrix)
  SE_signal_df[tmp_SE_name,] = tmp_column_means
}

    normalized_SE_signal_df_filter<- t(apply(SE_signal_df_filter, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
# tau
calculate_tau <- function(x) {
  n <- length(x)
  tau <- sum(1 - x) / (n - 1)
  return(tau)
}
tau_df  <- apply(normalized_SE_signal_df_filter, 1, calculate_tau) %>% as.data.frame();colnames(tau_df)[1] = 'tau_value'

con_1 = normalized_SE_signal_df_filter[,'SPP1+ Macro']>0.5
con_2 = tau_df$tau_value>0.7
normalized_SE_signal_df_filter_SPP1 = normalized_SE_signal_df_filter[con_1&con_2,]

for (i in 1:length(SE_gr)) {
  tmp_SE_gr = SE_gr[i]
  tmp_SE_name = paste0(c('SE',seqnames(tmp_SE_gr) %>% as.character(),ranges(tmp_SE_gr) %>% as.character()),collapse = '_')
  enhancer_hit = findOverlaps(fragment_info_gr,tmp_SE_gr) %>% queryHits(.)
  tmp_enhancer_name = fragment_info[enhancer_hit,'fragname']
  tmp_signal_matrix = signal_groupSE_mtx[tmp_enhancer_name,]
  tmp_column_means <- colMeans(tmp_signal_matrix)
  SE_signal_df[tmp_SE_name,] = tmp_column_means
}
  
}

{
  peak_mat <- centered_pmatPsB[region_vec,]
  # hierarchical cluster peaks and identify regulatory 'modules'
  hc <- hclust(dist(peak_mat), method="complete", members=NULL)
  clusters <- cutree(hc, h=23) # Cutoff of 25 yields 1-5 clusters per gene
  
  gene_pmat_list = list(mat=peak_mat[hc$order,], clusters=clusters[hc$order], hc=hc)

  hmat <- gene_pmat_list$mat
  clusters <- gene_pmat_list$clusters
  cmap <- getColorMap(cmaps_BOR$stallion, n=length(unique(clusters)))
  names(cmap) <- names(getFreqs(clusters))
  
  expr_df <- data.frame(
    geneExpr=L2gimatPsB[gene,],
    geneActivity=L2gsmatPsB[gene,],
    peakSum=log2(colSums(scaled_pmatPsB[rownames(hmat),])+1), # log2 of sum of scaled peaks
    cellType=psb_labels
  )
  
  peak_order <- rownames(hmat)
  exprRange <- quantile(expr_df$geneExpr, c(0, 0.95))
  use_cmap <- cmaps_BOR$sunrise
  col_fun <- circlize::colorRamp2(seq(exprRange[1], exprRange[2], length=length(use_cmap)), use_cmap)
  
  ht_opt$simple_anno_size <- unit(0.25, "cm")
  ta <- HeatmapAnnotation(cellType=expr_df$cellType, geneExpr=expr_df$geneExpr,
                          col=list(cellType=atacNamedClustCmap, geneExpr=col_fun))
  la <- HeatmapAnnotation(cluster=clusters[rownames(hmat)],col=list(cluster=cmap), 
                          which="row", show_legend=c("cluster"=TRUE))

hm <- BORHeatmap(
    hmat[peak_order,],
    limits=c(-2.,2.),
    clusterCols=TRUE, clusterRows=FALSE,
    labelCols=FALSE, labelRows=TRUE,
    dataColors = cmaps_BOR$solar_extra,
    top_annotation = ta,
    left_annotation = la,
    row_names_side = "left",
    border_gp=gpar(col="black"), 
    legendTitle="Row Z-score"
  )
  draw(hm)

  gene_pmat_list <- getLinkedPeakMat(gene, p2gGR, centered_pmatPsB)
  hmat <- gene_pmat_list$mat
  clusters <- gene_pmat_list$clusters
  cmap <- getColorMap(cmaps_BOR$stallion, n=length(unique(clusters)))
  names(cmap) <- names(getFreqs(clusters))
  

  expr_df <- data.frame(
    geneExpr=L2gimatPsB[gene,],
    geneActivity=L2gsmatPsB[gene,],
    peakSum=log2(colSums(scaled_pmatPsB[rownames(hmat),])+1), # 
    cellType=psb_labels
  )
  peak_rsq <- cor.test(expr_df$peakSum, expr_df$geneExpr)$estimate**2
  gs_rsq <- cor.test(expr_df$geneActivity, expr_df$geneExpr)$estimate**2
  
  yrange <- max(expr_df$geneExpr) - min(expr_df$geneExpr)
  plot_ylims <- c(min(expr_df$geneExpr) - 0.05*yrange, max(expr_df$geneExpr) + 0.1*yrange)
  xrange <- max(expr_df$peakSum) - min(expr_df$peakSum)
  plot_xlims <- c(min(expr_df$peakSum) - 0.05*xrange, max(expr_df$peakSum) + 0.1*xrange)
  
  pointSize <- 2.5
  
  expr_df_spp1 = expr_df %>% filter(cellType == 'SPP1+ Macro')
  

  fit <- lm(geneExpr ~ peakSum, data = expr_df_spp1)
  r2 <- summary(fit)$r.squared
  pval <- summary(fit)$coefficients[2, 4]

ggplot(expr_df_spp1, aes(x = peakSum, y = geneExpr)) +
    geom_point(aes(color = geneExpr), alpha = 1, size = 2) +
    scale_color_gradient(low = "darkorange", high = "#8B0000") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    annotate("text", x = max(expr_df$peakSum)*0.4, y = max(expr_df$geneExpr), 
             label = paste0("R² = ", round(r2, 2)), size = 4, hjust = 0) +
    annotate("text", x = max(expr_df$peakSum)*0.4, y = max(expr_df$geneExpr)*0.6, 
             label = paste0("P ", ifelse(pval < 0.001, "< 0.001", paste0("= ", signif(pval, 2)))), 
             size = 4, hjust = 0) +
    labs(
      x = "Log2(SE signal)",
      y = "Log2(SPP1 Expression)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    )
}

{
library(igraph)
library(tidygraph)
library(ggraph)
library(ggplot2)

calculateLinkageScore <- function(motifLocs, p2gGR){
    ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
    olGenes <- p2gGR[to(ol)]
    olGenes$motifScore <- motifLocs[from(ol)]$score
    olGenes$R2 <- olGenes$Correlation**2 
    LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
    LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
    LSdf$rank <- 1:nrow(LSdf)
    return(LSdf)
  }
  
  calculateMotifEnrichment <- function(motifLocs, p2gGR){
    motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
    m <- length(motifP2G) 
    n <- length(p2gGR) - m 
    motifLinks <- motifP2G$symbol %>% getFreqs()
    allLinks <- p2gGR$symbol %>% getFreqs()
    df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
    df$motifLinks[is.na(df$motifLinks)] <- 0
    df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
    df <- df[order(df$mLog10pval, decreasing=TRUE),]
    df$symbol <- rownames(df)
    return(df)
  }
  

get_circle_coords <- function(n, radius) {
  angles <- seq(0, 2 * pi, length.out = n + 1)[- (n + 1)]
  tibble(x = radius * cos(angles), y = radius * sin(angles))
}

layout_df <- bind_rows(
  get_circle_coords(length(TFs), 1),
  get_circle_coords(length(SEs), 2),
  get_circle_coords(length(Genes), 3)
)

ggraph(graph, layout = layout_df) +
  geom_edge_link(aes(edge_width = weight, edge_colour = edge_type, linetype = edge_type),
                 alpha = 0.6) +
  geom_node_point(aes(color = type), size = 3) +
  geom_node_text(aes(label = label), repel = TRUE, size = 2.5) +
  scale_color_manual(name = "Node Type",
                     values = c("TF" = "#FF6F61", "SE" = "#6BAED6", "Gene" = "#B2DF8A")) +
  scale_edge_colour_manual(name = "Edge Type",
                           values = c("TF_SE" = "#3182bd", "SE_Gene" = "gray40")) +
  scale_linetype_manual(name = "Edge Type",
                        values = c("TF_SE" = "dashed", "SE_Gene" = "solid")) +
  scale_edge_width(name = "Correlation", range = c(0.3, 1)) +
  theme_void() +
  theme(legend.position = "bottom")
}


{
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)

data <- read.delim("gene_interaction_results-2025_5_29 22_50_20.tsv", sep = "\t", header = TRUE)


data$interaction.score<- as.numeric(gsub("[^0-9\\.]", "", data$interaction.score))


gene_summary <- data %>%
  group_by(gene) %>%
  summarise(
    drug_count = n(),
    max_score = max(interaction.score, na.rm = TRUE)
  ) %>%
  filter(!is.na(max_score)) %>%     
  arrange(desc(max_score)) %>%
  mutate(rank = row_number())

gene_summary$max_score = log(gene_summary$max_score+1)

hr_data <- read_excel("HCC_cox_res.xlsx")

hr_data$gene_name <- toupper(hr_data$gene_name)

gene_summary_hr <- gene_summary %>%
  left_join(hr_data %>% select(gene = gene_name, HR = `hr...3`), by = "gene") %>%
  filter(!is.na(HR))  # 只保留有 HR 值的基因


pdf("drug_interact.pdf",height = 4.5,width = 6)
ggplot(gene_summary_hr, aes(x = rank, y = max_score)) +
  geom_point(aes(size = drug_count, fill = HR),
             shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  scale_fill_gradient(
    low = "white", high = "darkorange",
    limits = c(1, max(gene_summary_hr$HR, na.rm = TRUE))
  ) +
  scale_size(range = c(3, 12)) +
  geom_text_repel(
    aes(label = gene),
    size = 4,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    force = 2,
    force_pull = 0.5,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  labs(
    title = "Bubble Plot of Gene Interaction and Hazard Ratio",
    x = "Gene Rank by Max Interaction Score",
    y = "Max Interaction Score",
    fill = "Hazard Ratio (HR)",
    size = "Drug Count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black")
  )
dev.off()
}
