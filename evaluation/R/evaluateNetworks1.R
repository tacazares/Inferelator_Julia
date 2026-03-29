rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  library(circlize)
  # library(VennDiagram)
#   library(ggVennDiagram)
  library(grid)
  library(reshape2)
  library(RColorBrewer)
  library(purrr)

})

source("/data/miraldiNB/Michael/Scripts/VennDiagram.R")
source("/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/RprogUtils/evaluateNetUtils.R")

# dirOut <- "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/Bulk/5kbTSS/newPipe/spikeIgnored/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63"
dirOut <- "/data/miraldiNB/Michael/mCD4T_Wayman/Figures/Fig4/networkComparison1"
dir.create(dirOut, recursive = T, showWarnings = F)

# Combined/combined_max.tsv
# TFmRNA/edges_subset.txt
# netFiles <- list(
#                 "TFA" =    "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/TFA/edges_subset.tsv",
#                 "TFmRNA" = "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/TFmRNA/edges_subset.tsv",
#                 "Combined" = "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/Combined/combined_max.tsv"
#                 )

# netFiles <- list(
#                 "TFA" =    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.tsv",
#                 "TFmRNA" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.tsv",
#                 "Combined" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/combined_max.tsv"
#                 )

netFiles <- list(
            "PB" = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/combined_max.tsv",
            "SC" = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/noMergedTF/SC/ATAC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5/Combined/combined_max.tsv"
            # "MC2" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/metaCells/MC2/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63_logNorm/Combined/combined_max.tsv",
            # "SEA" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/metaCells/SEACells/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/combined_max.tsv"
)

# priorFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/MotifScan5kbTSS_b_sp.tsv"
priorFile = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10_sp.tsv"
# k <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/SCENICp/SCENICp_RE2Glinks_FIMO_5Kb_derived_b.tsv")
# k$Target <- rownames(k)
# b <- melt(k, id.vars = "Target", variable.names = "TF", value.name = "Weigths")  
# colnames(b)[2:3] <- c("TF", "Weights")   
# b <- b[, c("TF", "Target", "Weights")]  
# bb <- b[b$Weights != 0, ]  
# write.table(bb, "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/SCENICp/SCENICp_RE2Glinks_FIMO_5Kb_derived_sp.tsv", row.names = F, sep ="\t", quote = F)

tfCol <- "TF"        # specify TF column name here
targetCol <- "Gene"  # specify Target column name here
priorTfCol <- "TF"         # TF column name in prior file
priorTargetCol <- "Target" # Target column name in prior file
priorWgtCol <-"Weight"
nSelect <- 10
stabCol <- "Stability"
rankCol <- "signedQuantile"
compareNets <- TRUE
tfList <- NULL
lineageTFs <-c("Bcl6","Maf","Batf","Tox","Ascl2", "Foxp3","Ikzf2","Rorc", "Rorc","Stat3",
                 "Tbx21","Stat4","Eomes","Prdm1","Tcf7","Klf2","Lef1")

# lineageTFs <- c("TCF7","LEF1","ID3","KLF2","PRDM1","ZEB2","RUNX3","EOMES",
#                 "TBX21","STAT1","STAT4","HIF1A","RORC","RORA","STAT3","IRF4",
#                 "BATF","FOXP3","IKZF2","IKZF4","STAT5","FOXO1","GATA3","STAT6",
#                 "CIITA","RFX5","RFXAP","RFXANK","NLRC5")
N <- NULL
# file_k_clust_across <- "/Users/owop7y/Desktop/MiraldiLab/PROJECTS/GRN/Benchmark/mCD4T/Figures/Fig4/networkComparison/Kmeans_Jaccard_Edges_perTF.rds"
file_k_clust_across <- NULL

# ========================================================================
# ------- Read Input files
# ======================================================================== 

# ----- Read Prior
priorData <- read.table(priorFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(priorData, 3)
priorData <- priorData[priorData[[priorWgtCol]] != 0, ]
priorPairs <- unique(paste(priorData[[priorTfCol]], priorData[[priorTargetCol]], sep = "_"))

# --- Preload all network data into a list ---
netDataList <- list()
for (typeName in names(netFiles)) {
  filePath <- netFiles[[typeName]]
  netDataList[[typeName]] <- read.table(filePath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

numNetworks <- length(netDataList)
networkNames <- names(netDataList)

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ---- PART ONE:  Make a table of number of network ppts - uniqueTFs, uniqueTargets, Total Interactions, and % supported by prior
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# process each network
summaryList <- list()
for (typeName in networkNames) {
  data <- netDataList[[typeName]]

  uniqueTf <- length(unique(data[[tfCol]]))
  uniqueTarget <- length(unique(data[[targetCol]]))

  pairStrings <- unique(paste(data[[tfCol]], data[[targetCol]], sep = "_"))
  totalInteractions <- length(pairStrings)
  
  supportedCount <- length(intersect(pairStrings, priorPairs))
  percentSupportedByPrior <- 100 * supportedCount / totalInteractions

  summaryList[[typeName]] <- data.frame(type = typeName, uniqueTf = uniqueTf, 
                            uniqueTarget = uniqueTarget, totalInteractions = totalInteractions,
                            percentSupportedByPrior = percentSupportedByPrior)
}   

summaryTable <- do.call(rbind, summaryList)
print(summaryTable)

write.table(as.data.frame(summaryTable), file.path(dirOut, "summaryStatistics.tsv"), col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ---- PART TWO:  Compare Networks (can handle more than 2)
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

# ──── Compare TF-Targets across methods
# Initialize a data frame to store results
tfTargetCounts <- data.frame(TF = unique(unlist(lapply(netDataList, function(x) unique(x$TF)))))

# Count targets for each TF in each network
for(netName in names(netDataList)) {
  netData <- netDataList[[netName]]
  
  # Count number of unique target genes per TF
  targetCounts <- netData %>%
    group_by(TF) %>%
    summarise(targetCount = n())
  
  colname <- netName
  tfTargetCounts <- tfTargetCounts %>%
    left_join(targetCounts, by = "TF") %>%
    rename(!!colname := targetCount)
}
# Replace NA with 0 for TFs not present in a network
tfTargetCounts[is.na(tfTargetCounts)] <- 0

# Build file suffix
fileSuffix <- if(!is.null(topNpercent)) {
    paste0(rankCol, "_top", topNpercent, "pct")
} else {
    paste0(rankCol, "_all")
}
networkNames <- names(netDataList)
numNetworks <- length(netDataList)

edgeSets <- getEdgeSets(netDataList, tfCol, targetCol, rankCol)

# --------  unique-edge counts ----------------------------------------- ##
uniqCountsDf <- data.frame(
        network = names(edgeSets),
        numUniqueEdges = sapply(seq_along(edgeSets), function(k){
                                    curr <- edgeSets[[k]]
                                    others <- unique(unlist(edgeSets[-k]))
                                    sum(!curr %in% others)
                                }),
        row.names = NULL,
        stringsAsFactors = FALSE
    )

write.table(uniqCountsDf,
                        file.path(dirOut, paste0("num_UniqueEdges_", fileSuffix, ".tsv")),
                        quote=F, col.names=NA, row.names=TRUE, sep="\t")

# -------  Plot global Jaccard heatmap
# pairs <- plotGlobalJaccard(edgeSets, fileOut = file.path(dirOut, paste0("Global_Jaccard_",fileSuffix, ".pdf")), fig_width = 2.5, fig_height = 2.5)
# pairInt <- pairs[["pairInt"]]

plotSimilarityHeatmap(pairs[["pairJac"]], fileOut=file.path(dirOut, paste0("Global_Jaccard_",fileSuffix, ".pdf")), fontsize_number=7, fontsize=9, fig_width=2.5, fig_height=2.5) 
pairDf <- melt(pairs[["pairInt"]], varnames = c("network1", "network2"), value.name = "numIntersection", na.rm = TRUE)  
write.table(pairDf,file.path(dirOut, paste0("num_Pairwise_Intersection_", fileSuffix, ".tsv")), quote=F, col.names=NA, row.names=TRUE, sep="\t")


# --------- Compute Jaccard for multiple topN% -------------
topNpercentList <- seq(10, 100, by = 10)
jaccardResults <- list()
for (topN in topNpercentList) {
    currEdgeSets <- getEdgeSets(netDataList, tfCol, targetCol, rankCol, "topNpercent", topN)
    # Compute pairwise Jaccard
    fileSuffix <- paste0(rankCol, "_top", topN, "pct")
    pairs <- computeJaccardMatrix(currEdgeSets)
    pairJac <- pairs[["pairJac"]]
    # Convert pairwise matrix to long format
    pairJacDf <- as.data.frame(as.table(pairJac))
    colnames(pairJacDf) <- c("Network1", "Network2", "Jaccard")
    pairJacDf$TopNpercent <- topN
    jaccardResults[[as.character(topN)]] <- pairJacDf
}
# Combine all percentages
jaccardDf <- do.call(rbind, jaccardResults)
# Filter out self-comparisons
jaccardDf <- jaccardDf[jaccardDf$Network1 != jaccardDf$Network2, ]
jaccardDf<- jaccardDf %>%
    dplyr::mutate(
      Network1 = as.character(Network1),
      Network2 = as.character(Network2),
      pairID = ifelse(Network1 < Network2,
                      paste(Network1, Network2, sep = "-"),
                      paste(Network2, Network1, sep = "-"))
    ) %>%
    dplyr::group_by(TopNpercent, pairID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
write.table(jaccardDf,file.path(dirOut, "Jaccard_Sim.tsv"), quote=F, col.names=NA, row.names=TRUE, sep="\t")

# Create color mapping based on Set1 palette
uniquePairs <- sort(unique(jaccardDf$pairID))
colors <- c("#449B75" ,"#FFE528",  "#999999", "#E41A1C", "#AC5782", "#C66764")
comparisonColor = setNames(colors, uniquePairs)

# Plot line graph with ggplot2
p_jac <- ggplot(jaccardDf, aes(x = TopNpercent, y = Jaccard,color = pairID,group = pairID)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    scale_color_manual(values = comparisonColor) +
    scale_x_continuous(breaks = seq(10, 100, by = 10)) +
    scale_y_continuous(limits = c(0,1)) +
    labs(
      x = "Top N% edges",
      y = "Jaccard similarity",
      color = "Network pair"
    ) +
    theme_bw(base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.box.just = "left",   # centers the legend
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.1)
  )
ggsave(file.path(dirOut, "linePlot_JaccardSim.pdf"),plot = p_jac, width = 4, height = 3, dpi = 600)

#plotEdgeSharing(edgeSets, fileSuffix, outDir) 
                                        
# --- 3. TF-centric analysis
# For each TF, how consistent are its predicted targets across networks?
# Biological networks are TF-driven, so checking target consistency per TF
tfMatrix <- computeTFJaccard(netDataList, tfList)
# MAke Heatmap
getPalette = colorRampPalette(brewer.pal(13, "Set1"))
comparisonColor = getPalette(length(colnames(tfMatrix)))
comparisonColor = setNames(comparisonColor, colnames(tfMatrix))
colAnn <- HeatmapAnnotation(
                `Comparison` = colnames(tfMatrix),
                col = list('Comparison' = comparisonColor),
                annotation_legend_param = list(
                    Comparison = list(
                        title_gp = gpar(fontsize = 8, fontface = "plain"),   # title size
                        labels_gp = gpar(fontsize = 6)                       # category label size
                    )),    
                    simple_anno_size = unit(3, "mm"),
                      show_annotation_name = FALSE
                    )   


tfMatrix[is.na(tfMatrix)] <- 0
tfRobustness <- rowMeans(tfMatrix, na.rm = TRUE)  # average Jaccard per TF

k_center = 6
if (is.null(file_k_clust_across)){
  k_clust_across <- kmeans(tfMatrix, centers=k_center,nstart=20, iter.max = 50)
}else{
  k_clust_across <- readRDS(file_k_clust_across)
}
split_by <- factor(k_clust_across$cluster, levels = c(1:6))

# Add row markers (labels with lines)
row_mark <- rowAnnotation(
  mark = anno_mark(
    at = which(rownames(tfMatrix) %in% lineageTFs),
    labels = lineageTFs,
    labels_gp = gpar(fontsize = 6, fontface = "plain"),
    padding = unit(0.3, "mm"),
    side = "left",        # line starts from left of heatmap
  )
)

annot_cols <- c('1'='#fb940b','2'='#01cc00','3'='#2085ec','4'='#fe98bf','5'='#ad7a5b', '6'='turquoise4')
df_clust <- data.frame(across=k_clust_across$cluster)
df_clust$across <- factor(df_clust$across, levels=1:k_center, ordered=T)
rowAnn_left <-  HeatmapAnnotation(`Cluster`= df_clust$across,
                                    col = list('Cluster'= annot_cols),
                                    which = 'row',
                                    simple_anno_size = unit(2, 'mm'),
                                    show_legend = FALSE,
                                    show_annotation_name = F)

# robustCols <- colorRamp2(c(min(tfRobustness), max(tfRobustness)), c("lightgrey","blue"))
# robustCols <- colorRamp2( seq(from = min(tfRobustness), to = max(tfRobustness), length.out = 200),inferno(200))
robustCols <- colorRamp2(seq(min(tfRobustness), max(tfRobustness), length.out = 100), colorRampPalette(brewer.pal(9, "Blues"))(100))
rowRobustness <- rowAnnotation(
  Robustness = anno_simple(
    tfRobustness,
    col = robustCols,,
    border = TRUE
  ),
  show_annotation_name = F,
  width = unit(1.5, "mm")
)

# Create a separate legend for robustness
robustLegend <- Legend(
  title = "TF Robustness",
  col_fun = robustCols,
  title_gp = gpar(fontsize = 8, fontface = "plain"),   # unbold
  labels_gp = gpar(fontsize = 6)                        # optional: label size
)

# dirOut <- "/data/MiraldiLab/team/Michael/mCD4T_WaymanFigures/Fig4/networkComparison/"
# dir.create(dirOut, recursive = T, showWarnings = F)

heatCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(100)
pdf(file.path(dirOut,  "Jaccard_Edges_perTF1.pdf"), width = 2.3, height = 4.5, compress = T)
ht1 <- Heatmap(tfMatrix,
            name='Jaccard Similarity',
            show_row_names = F,
            row_split=split_by,
            row_gap = unit(0, "mm"),
            column_gap = unit(0, "mm"), 
            border = TRUE,
            row_title = NULL,
            row_dend_reorder = F,
            use_raster = F,
            col = heatCols,
            # clustering settings
            cluster_rows = TRUE,                    # allow hierarchical clustering
            cluster_row_slices = TRUE,              # reorder within each k-means cluster
            cluster_columns = FALSE,
            cluster_column_slices = F,
            # column_order = colnames(tfa_z_across),
            show_row_dend = FALSE,                   # show dendrogram within slices
            show_column_dend = FALSE,
            show_column_names = FALSE,
            column_names_side = 'top',
            # column_names_rot = 45,
            column_title = NULL,
            column_split = colnames(tfMatrix),
            top_annotation = colAnn,
            left_annotation = rowAnn_left,
            heatmap_legend_param = list(
                    direction = "horizontal",
                    title = "Jaccard",
                    title_position = "topcenter",
                    title_gp = gpar(fontsize = 8, fontface = "plain"),   # legend title size
                    labels_gp = gpar(fontsize = 6),                     # numbers/tick labels size
                    legend_width = unit(2.5, "cm"),
                    legend_height = unit(0.5, "cm")
                )
            )   
draw(row_mark+ht1+rowRobustness , heatmap_legend_side = "bottom")
dev.off()

ht2 <- draw(ht1)
row_order_list <- row_order(ht2)
saveRDS(row_order_list, file.path(dirOut, "final_row_order.rds") )
# Save Legend as PDF
pdf(file.path(dirOut, "TF_Robustness_Legend.pdf"), width = 2, height = 2)
draw(robustLegend)
dev.off()
saveRDS(k_clust_across, file.path(dirOut, "Kmeans_Jaccard_Edges_perTF.rds"))


# ---------- Check 
clust_df <- df_clust
clust_df$TF <- rownames(clust_df)
colnames(clust_df)[1] <- "cluster"
tfTargetCounts <- tfTargetCounts %>%
  left_join(clust_df, by = "TF")

# Calculate average targets across all methods
tfTargetCounts$avgTargets <- rowMeans(tfTargetCounts[, names(netDataList)])
tfTargetCounts$medianTargets <- apply(tfTargetCounts[, names(netDataList)],1,median, na.rm = T)
# # Summary statistics by cluster
# clusterSummary <- tfTargetCounts %>%
#   group_by(cluster) %>%
#   summarise(
#     meanTargets = mean(avgTargets),
#     medianTargets = median(avgTargets),
#     sdTargets = sd(avgTargets),
#     nTFs = n()
#   )

write.table(tfTargetCounts, file.path(dirOut, "TF_avgTargetCounts_byRobustnessCluster.tsv"), quote=F, col.names=NA, row.names=TRUE, sep="\t")

pAvg <- ggplot(tfTargetCounts, aes(x = cluster, y = avgTargets, fill = cluster)) +
  geom_boxplot(outlier.color = "red", outlier.size = 1,outlier.shape = 19, outlier.alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.1) +
  labs(title = "",
       y = "Average # of Targets",
       x = "") +
  scale_fill_manual(values = annot_cols) +
  theme_minimal() +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank() ,
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1)
  )
ggsave(file.path(dirOut, "TF_avgTargetCounts_Distribution.pdf"), pAvg, width = 2.5, height = 2, dpi = 600)


pMedian <- ggplot(tfTargetCounts, aes(x = cluster, y = medianTargets, fill = cluster)) +
  geom_boxplot(outlier.color = "red", outlier.size = 1,outlier.shape = 19, outlier.alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.1) +
  labs(title = "",
       y = "Median # of Targets",
       x = "") +
  scale_fill_manual(values = annot_cols) +
  theme_minimal() +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank() ,
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1)
  )
ggsave(file.path(dirOut, "TF_medianTargetCounts_Distribution.pdf"), pMedian, width = 2.5, height = 2, dpi = 600)


# tf_stats <- lapply(names(netDataList), function(nm){
#   df <- netDataList[[nm]]
#   df %>%
#     count(!!sym(tfCol)) %>%
#     rename(TF = !!sym(tfCol), !!nm := n)
# }) %>%
#   reduce(full_join, by = "TF") %>%
#   mutate(cluster = k_clust_across$cluster[TF])
# 
# # Which cluster have higher mean degree and which has lower mean degree
# cluster_summary <- tf_stats %>%
#   group_by(cluster) %>%
#   summarise(across(starts_with("PB") | starts_with("SC") | starts_with("MC") | starts_with("SEA"),
#                    list(mean = ~mean(.x, na.rm=TRUE),
#                         sd = ~sd(.x, na.rm=TRUE))),
#             n_TFs = n())
# 
# ## Only Lineage TFs
# tfMatrix_lineage <- tfMatrix[intersect(lineageTFs, rownames(tfMatrix)), ]
# k_clust_lineage<- kmeans(tfMatrix_lineage, centers=4,nstart=20, iter.max = 50)
# split_by_lineage <- factor(k_clust_lineage$cluster, levels = c(1:4))
# 
# pdf(file.path(dirOut,  "Jaccard_Edges_perTF_lineage.pdf"), width = 2, height = 3, compress = T)
# ht <- Heatmap(tfMatrix_lineage,
#                name = 'Jaccard Similarity',
#                show_row_names = TRUE,
#                row_names_gp = gpar(fontsize = 6),
#                row_names_side = "left",           # rownames on the left
#                row_split = split_by_lineage,
#                row_gap = unit(0, "mm"),
#                column_gap = unit(0, "mm"),
#                border = TRUE,
#                row_title = NULL,
#                row_dend_reorder = FALSE,
#                use_raster = FALSE,
#                col = heatCols,
#                cluster_rows = TRUE,
#                cluster_row_slices = TRUE,
#                cluster_columns = FALSE,
#                cluster_column_slices = FALSE,
#                show_row_dend = FALSE,
#                show_column_dend = FALSE,
#                show_column_names = FALSE,
#                column_names_side = 'top',
#                column_title = NULL,
#                column_split = colnames(tfMatrix_lineage),
#                top_annotation = colAnn,
#                heatmap_legend_param = list(
#                     direction = "horizontal",
#                     title = "Jaccard Similarity",
#                     title_position = "topcenter",
#                     title_gp = gpar(fontsize = 8, fontface = "plain"),   # legend title size
#                     labels_gp = gpar(fontsize = 6),                     # numbers/tick labels size
#                     legend_width = unit(2.5, "cm"),
#                     legend_height = unit(0.5, "cm")
#                 )
# )
# draw(ht, heatmap_legend_side = "bottom")
# dev.off()
# saveRDS(k_clust_lineage, file.path(dirOut, "Kmeans_Jaccard_Edges_perTF_lineage.rds"))

# --- 4. Aggregated network-pair Jaccard
fun = median
fun_name = "median"  # returns "mean"
aggMat <- computeAggregatedPairJaccard(tfMatrix, fun = fun)
plotAggregatedJaccard(aggMat,file.path(dirOut, paste0("AggregatedTF_Jaccard_", fun_name, ".pdf")), fig_width = 2.5, fig_height = 2.5)

fun = mean
fun_name = "mean"  # returns "mean"
aggMat <- computeAggregatedPairJaccard(tfMatrix, fun = fun)
plotAggregatedJaccard(aggMat,file.path(dirOut, paste0("AggregatedTF_Jaccard_", fun_name, ".pdf")), fig_width = 2.5, fig_height = 2.5)


# ----- Are top regulators (high-degree TFs) stable across pseudobulk / single-cell / metacell networks.
hubSim <- computeHubOverlapHeatmap(
  netDataList = netDataList, tfCol = tfCol, networkNames = networkNames,
  dirOut = dirOut, topN = 50,
  metric = "jaccard",   # or "overlap"
  fontsize = 9, fontsize_number = 7, fig_width = 2.5, fig_height = 2.5
)

hubSim <- computeHubOverlapHeatmap(
  netDataList = netDataList, tfCol = tfCol, networkNames = networkNames,
  dirOut = dirOut, topN = 50,
  metric = "overlap",   # or "jaccard"
  fontsize = 9, fontsize_number = 7, fig_width = 2.5, fig_height = 2.5
)

# ────────────────────────────────────────────────────────────────────────────────────────────
# PART THREE:   Histogram distributions of:
                # A: Targets per TF: # times each gene is targeted (number of TFs regulating it).
                # B: TFs per Target: # times each TF is a regulator (number of genes it controls)
                # C: Box-Plot of Stability of the top N low and high in degree genes

                #NOTE:  Each Figure is saved in the same path as the network being evaluated
# ────────────────────────────────────────────────────────────────────────────────────────────
# Generate both plots for each network
for (typeName in names(netDataList)) {
  data <- netDataList[[typeName]]
  basePath <- dirname(netFiles[[typeName]])

  # Plot 1: # Targets per TF 
  # tfTargetCounts <- table(data[[tfCol]])
  # dfTF <- data.frame(TF = names(tfTargetCounts), targetCount = as.integer(tfTargetCounts))
  dfTF <- tfTargetCounts[, c("TF", tfCol), drop = F]
  colnames(dfTF)[2] <- "targetCount"

  p1 <- ggplot(dfTF, aes(x = targetCount)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "#0072B2", color = "black", alpha = 0.8) +
    # scale_x_continuous(breaks = seq(0, max(dfTF$targetCount), 1)) +
    labs(title = paste("Distribution of # Targets per TF -", typeName),
         x = "# Targets", y = "# TFs") +
    theme_bw(base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # axis.line = element_line(color = "black", linewidth = 0.4),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5),
      # panel.background = element_rect("black", fill = NA)
    )

  ggsave(file.path(basePath, paste0("TargetCountPerTF_", typeName, ".pdf")),
         plot = p1, width = 3, height = 3)

  # Plot 2: # TFs per Target
  targetTFCounts <- table(data[[targetCol]])
  dfTarget<- data.frame(Target = names(targetTFCounts), tfCount = as.integer(targetTFCounts))
  dfTargetSorted <- dfTarget[order(dfTarget$tfCount), ]

  p2 <- ggplot(dfTarget, aes(x = tfCount)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    # scale_x_continuous(breaks = seq(0, max(dfTarget$tfCount), 1)) +
    labs(title = paste("Distribution of # TFs per Target -", typeName),
         x = "# TFs", y = "# Targets") +
    theme_bw(base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # axis.line = element_line(color = "black", linewidth = 0.4),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5),
      # panel.background = element_rect("black", fill = NA)
    )

  ggsave(file.path(basePath, paste0("TFCountPerTarget_", typeName, ".pdf")),
         plot = p2, width = 3, height = 3)
  
  # Select top N low and high TF targets
  lowTargs <- dfTargetSorted$Target[1:nSelect]
  highTargs <- tail(dfTargetSorted, nSelect)$Target

  stabilityDF <- data[data[[targetCol]] %in% c(lowTargs, highTargs),  c(targetCol, stabCol)]
  # Add Group column based on whether the target is in lowTargs or highTargs
 stabilityDF$Group <- ifelse(stabilityDF[[targetCol]] %in% lowTargs,
                            "Low TFs per Target",
                            "High TFs per Target")
  stabilityDF <- merge(stabilityDF, dfTargetSorted, by.x = targetCol, by.y = "Target")
  stabilityDF$targetLabel <- paste0(stabilityDF[[targetCol]], "(", stabilityDF$tfCount, ")")
  # Keep plotting order
#   stabilityDF$targetLabel <- factor(stabilityDF$targetLabel,
#                                      levels = unique(stabilityDF$targetLabel))
                                     # Order based on tfCount directly
    stabilityDF$targetLabel <- factor(
    stabilityDF$targetLabel,
    levels = unique(stabilityDF$targetLabel[order(stabilityDF$tfCount)])
    )

  # Plot: boxplot per target
  p3 <- ggplot(stabilityDF, aes(x = targetLabel, y = Stability, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    labs(title = paste("Per-Target Stability Distribution -", typeName),
         x = "Number of TFs per Target (TF count)", y = "Stability") +
    theme_minimal() +
    theme_bw(base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # axis.line = element_line(color = "black", linewidth = 0.4),
      axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "top"
      # panel.background = element_rect("black", fill = NA)
    )
  
  # Save
  ggsave(file.path(basePath, paste0("Top", nSelect, "HighorLow_inDegreeGenes_Boxplot", typeName, ".pdf")),
         plot = p3, width = 12, height = 5)

}




# 



# 
# library(reshape2)
# library(ggplot2)
# 
# # Convert to long format for ggplot
# df_long <- melt(tfa_by_celltype)
# colnames(df_long) <- c("TF", "CellType", "TFA")
# 
# # Add a flag for lineage TFs
# df_long$Lineage <- mapply(function(tf, ct) {
#   if(tf %in% lineage_TFs[[ct]]) "Lineage" else "Other"
# }, df_long$TF, df_long$CellType)
# 
# # Boxplot or violin plot
# ggplot(df_long, aes(x = CellType, y = TFA, fill = Lineage)) +
#   geom_violin(alpha = 0.6) +
#   geom_boxplot(width=0.1, outlier.shape=NA) +
#   scale_fill_manual(values = c("Lineage"="red", "Other"="grey")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   labs(title="Lineage TF enrichment in TFA", y="Mean TFA")
# 
# 
# lineage_TFs <- list(
#   Tfh10   = c("Bcl6","Maf","Batf"),
#   Tfh_Int = c("Bcl6","Maf","Batf"),
#   Tfh     = c("Bcl6","Maf"),
#   Tfr     = c("Foxp3","Bcl6"),
#   cTreg   = c("Foxp3","Ikzf2"),
#   eTreg   = c("Foxp3","Ikzf2"),
#   rTreg   = c("Foxp3","Ikzf2"),
#   Treg_Rorc = c("Foxp3","Rorc"),
#   Th17    = c("Rorc","Batf","Stat3"),
#   Th1     = c("Tbx21","Stat4"),
#   CTL_Prdm1 = c("Tbx21","Eomes","Prdm1"),
#   CTL_Bcl6 = c("Tbx21","Eomes","Prdm1"),
#   TEM     = c("Eomes","Tcf7","Klf2"),
#   TCM     = c("Eomes","Tcf7","Klf2")
# )
# 

# lineage_tfs_list <- list(
#   TCM = c("TCF7", "LEF1", "ID3", "KLF2"),
#   TEM = c("PRDM1", "ZEB2", "RUNX3", "EOMES"),
#   Th1 = c("TBX21", "STAT1", "STAT4", "RUNX3", "HIF1A"),
#   Th17 = c("RORC", "RORA", "STAT3", "IRF4", "BATF"),
#   Treg = c("FOXP3", "IKZF2", "IKZF4", "STAT5"),
#   Naive = c("TCF7", "LEF1", "KLF2", "FOXO1"),
#   Th2 = c("GATA3", "STAT5", "STAT6", "IRF4"),
#   CTL = c("EOMES", "TBX21", "RUNX3", "ZEB2", "PRDM1"),
#   MHCII = c("CIITA", "RFX5", "RFXAP", "RFXANK", "NLRC5")
# )

# 
# 
# #  Distribution per F1.
# # Also show F1 score based on heatmap
# ggplot(pr_df, aes(x = network, y = F1)) +
#   geom_boxplot(outlier.size = 0.8) +
#   geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1)) +
#   ylab("F1 score") + xlab("Representation")
