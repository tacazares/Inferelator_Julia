rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  library(circlize)
  library(grid)
  library(reshape2)
  library(RColorBrewer)
  library(purrr)
})

source("/path/to/evaluation/R/evaluateNetUtils.R")

# ========================================================================
# USER INPUTS — edit this section before running
# ========================================================================

# Output directory
dirOut <- "/path/to/output/networkComparison"
dir.create(dirOut, recursive = TRUE, showWarnings = FALSE)

# Named list of networks to compare: label => path to edges file
# Typical files: edges.tsv, edges_subset.tsv, combined_max.tsv
netFiles <- list(
  "NetworkA" = "/path/to/networkA/Combined/combined_max.tsv",
  "NetworkB" = "/path/to/networkB/Combined/combined_max.tsv"
)

# Prior network file (sparse format: TF, Target, Weight columns)
priorFile <- "/path/to/priors/prior_sp.tsv"

# Column names in edge files
tfCol     <- "TF"              # TF column
targetCol <- "Gene"            # Target gene column
rankCol   <- "signedQuantile"  # Ranking score column
stabCol   <- "Stability"       # Stability column

# Prior column names
priorTfCol     <- "TF"
priorTargetCol <- "Target"
priorWgtCol    <- "Weight"

# Analysis parameters
nSelect            <- 10    # Number of top/bottom degree genes for stability boxplot
N                  <- NULL  # Top-N% edges to use for comparisons (NULL = all edges)
tfList             <- NULL  # Optional: restrict TF Jaccard to a subset of TFs (NULL = all)
k_center           <- 6     # Number of k-means clusters for TF Jaccard heatmap
file_k_clust       <- NULL  # Path to saved k-means RDS to reuse clustering (NULL = recompute)

# Optional: marker TFs to annotate in the TF Jaccard heatmap (set to c() to skip)
lineageTFs <- c("TF1", "TF2", "TF3")

# ========================================================================
# LOAD DATA
# ========================================================================

priorData  <- read.table(priorFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
priorData  <- priorData[priorData[[priorWgtCol]] != 0, ]
priorPairs <- unique(paste(priorData[[priorTfCol]], priorData[[priorTargetCol]], sep = "_"))

netDataList  <- lapply(netFiles, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
networkNames <- names(netDataList)
numNetworks  <- length(netDataList)

# ========================================================================
# PART 1: Summary statistics
# How many TFs, targets, and interactions does each network have?
# How many edges are supported by the prior?
# Output: summaryStatistics.tsv
# ========================================================================

summaryList <- lapply(networkNames, function(nm) {
  data        <- netDataList[[nm]]
  pairStrings <- unique(paste(data[[tfCol]], data[[targetCol]], sep = "_"))
  supported   <- length(intersect(pairStrings, priorPairs))
  data.frame(
    Network             = nm,
    uniqueTFs           = length(unique(data[[tfCol]])),
    uniqueTargets       = length(unique(data[[targetCol]])),
    totalInteractions   = length(pairStrings),
    pctSupportedByPrior = round(100 * supported / length(pairStrings), 2)
  )
})
summaryTable <- do.call(rbind, summaryList)
print(summaryTable)
write.table(summaryTable, file.path(dirOut, "summaryStatistics.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# ========================================================================
# PART 2: Pairwise network comparison
# How similar are the networks globally and across confidence thresholds?
# Outputs: Global_Jaccard_*.pdf, num_Pairwise_Intersection_*.tsv,
#          num_UniqueEdges_*.tsv, Jaccard_Sim.tsv, linePlot_JaccardSim.pdf
# ========================================================================

# Count targets per TF across all networks (used later in Part 4)
tfTargetCounts <- data.frame(TF = unique(unlist(lapply(netDataList, function(x) unique(x[[tfCol]])))))
for (nm in networkNames) {
  counts <- netDataList[[nm]] %>%
    group_by(.data[[tfCol]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    rename(TF = 1)
  tfTargetCounts <- tfTargetCounts %>%
    left_join(counts, by = "TF") %>%
    rename(!!nm := n)
}
tfTargetCounts[is.na(tfTargetCounts)] <- 0

edgeSets   <- getEdgeSets(netDataList, tfCol, targetCol, rankCol, N = N)
fileSuffix <- if (!is.null(N)) paste0(rankCol, "_top", N, "pct") else paste0(rankCol, "_all")

# Unique edges per network
uniqCountsDf <- data.frame(
  network        = names(edgeSets),
  numUniqueEdges = sapply(seq_along(edgeSets), function(k) {
    sum(!edgeSets[[k]] %in% unique(unlist(edgeSets[-k])))
  }),
  stringsAsFactors = FALSE
)
write.table(uniqCountsDf,
            file.path(dirOut, paste0("num_UniqueEdges_", fileSuffix, ".tsv")),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

# Global Jaccard heatmap
pairs  <- computeJaccardMatrix(edgeSets)
plotSimilarityHeatmap(pairs[["pairJac"]],
                      fileOut    = file.path(dirOut, paste0("Global_Jaccard_", fileSuffix, ".pdf")),
                      fig_width  = 2.5, fig_height = 2.5)
pairDf <- melt(pairs[["pairInt"]], varnames = c("network1", "network2"),
               value.name = "numIntersection", na.rm = TRUE)
write.table(pairDf,
            file.path(dirOut, paste0("num_Pairwise_Intersection_", fileSuffix, ".tsv")),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

# Jaccard across top-N% thresholds
topNpercentList <- seq(10, 100, by = 10)
jaccardResults  <- lapply(topNpercentList, function(topN) {
  currEdgeSets <- getEdgeSets(netDataList, tfCol, targetCol, rankCol, "topNpercent", topN)
  pairJac      <- computeJaccardMatrix(currEdgeSets)[["pairJac"]]
  df           <- as.data.frame(as.table(pairJac))
  colnames(df) <- c("Network1", "Network2", "Jaccard")
  df$TopNpercent <- topN
  df
})
jaccardDf <- do.call(rbind, jaccardResults)
jaccardDf <- jaccardDf[jaccardDf$Network1 != jaccardDf$Network2, ] %>%
  mutate(pairID = ifelse(Network1 < Network2,
                         paste(Network1, Network2, sep = "-"),
                         paste(Network2, Network1, sep = "-"))) %>%
  group_by(TopNpercent, pairID) %>% slice(1) %>% ungroup()
write.table(jaccardDf, file.path(dirOut, "Jaccard_Sim.tsv"),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

uniquePairs     <- sort(unique(jaccardDf$pairID))
comparisonColor <- setNames(
  colorRampPalette(brewer.pal(min(length(uniquePairs) + 1, 9), "Set1"))(length(uniquePairs)),
  uniquePairs
)
p_jac <- ggplot(jaccardDf, aes(x = TopNpercent, y = Jaccard, color = pairID, group = pairID)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1) +
  scale_color_manual(values = comparisonColor) +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Top N% edges", y = "Jaccard similarity") +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text  = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.1))
ggsave(file.path(dirOut, "linePlot_JaccardSim.pdf"), plot = p_jac,
       width = 4, height = 3, dpi = 600)

# ========================================================================
# PART 3: TF-centric Jaccard
# For each TF, how consistent are its predicted targets across networks?
# High Jaccard = robust TF regulon; Low Jaccard = network-specific.
# Outputs: Jaccard_Edges_perTF.pdf, AggregatedTF_Jaccard_*.pdf,
#          HubTF_Overlap_*.pdf, TF_avgTargetCounts_*.pdf,
#          Kmeans_Jaccard_Edges_perTF.rds
# ========================================================================

tfMatrix     <- computeTFJaccard(netDataList, tfList)
tfMatrix[is.na(tfMatrix)] <- 0
tfRobustness <- rowMeans(tfMatrix, na.rm = TRUE)

# K-means clustering (reuse saved clustering if provided)
k_clust  <- if (is.null(file_k_clust)) {
  kmeans(tfMatrix, centers = k_center, nstart = 20, iter.max = 50)
} else {
  readRDS(file_k_clust)
}
split_by <- factor(k_clust$cluster, levels = seq_len(k_center))

# Color setup
heatCols    <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
robustCols  <- colorRamp2(seq(min(tfRobustness), max(tfRobustness), length.out = 100),
                          colorRampPalette(brewer.pal(9, "Blues"))(100))
nPairs      <- ncol(tfMatrix)
compColor   <- setNames(
  colorRampPalette(brewer.pal(min(nPairs + 1, 13), "Set1"))(nPairs),
  colnames(tfMatrix)
)
annotCols <- setNames(
  brewer.pal(max(k_center, 3), "Set2")[seq_len(k_center)],
  as.character(seq_len(k_center))
)

colAnn <- HeatmapAnnotation(
  Comparison = colnames(tfMatrix),
  col = list(Comparison = compColor),
  annotation_legend_param = list(
    Comparison = list(title_gp  = gpar(fontsize = 8, fontface = "plain"),
                      labels_gp = gpar(fontsize = 6))),
  simple_anno_size = unit(3, "mm"), show_annotation_name = FALSE
)
df_clust    <- data.frame(Cluster = factor(k_clust$cluster, levels = seq_len(k_center)))
rowAnn_left <- HeatmapAnnotation(
  Cluster = df_clust$Cluster,
  col = list(Cluster = annotCols), which = "row",
  simple_anno_size = unit(2, "mm"), show_legend = FALSE, show_annotation_name = FALSE
)
rowRobustness <- rowAnnotation(
  Robustness = anno_simple(tfRobustness, col = robustCols, border = TRUE),
  show_annotation_name = FALSE, width = unit(1.5, "mm")
)
row_mark <- if (length(lineageTFs) > 0 && any(lineageTFs %in% rownames(tfMatrix))) {
  rowAnnotation(mark = anno_mark(
    at     = which(rownames(tfMatrix) %in% lineageTFs),
    labels = rownames(tfMatrix)[rownames(tfMatrix) %in% lineageTFs],
    labels_gp = gpar(fontsize = 6, fontface = "plain"),
    padding = unit(0.3, "mm"), side = "left"
  ))
} else { NULL }

robustLegend <- Legend(title = "TF Robustness", col_fun = robustCols,
                       title_gp  = gpar(fontsize = 8, fontface = "plain"),
                       labels_gp = gpar(fontsize = 6))

pdf(file.path(dirOut, "Jaccard_Edges_perTF.pdf"), width = 2.3, height = 4.5, compress = TRUE)
ht <- Heatmap(tfMatrix, name = "Jaccard",
              show_row_names = FALSE,
              row_split = split_by, row_gap = unit(0, "mm"),
              column_gap = unit(0, "mm"), border = TRUE,
              row_title = NULL, row_dend_reorder = FALSE, use_raster = FALSE,
              col = heatCols,
              cluster_rows = TRUE, cluster_row_slices = TRUE,
              cluster_columns = FALSE, cluster_column_slices = FALSE,
              show_row_dend = FALSE, show_column_dend = FALSE,
              show_column_names = FALSE, column_names_side = "top",
              column_title = NULL, column_split = colnames(tfMatrix),
              top_annotation = colAnn, left_annotation = rowAnn_left,
              heatmap_legend_param = list(
                direction = "horizontal", title = "Jaccard",
                title_position = "topcenter",
                title_gp  = gpar(fontsize = 8, fontface = "plain"),
                labels_gp = gpar(fontsize = 6),
                legend_width = unit(2.5, "cm"), legend_height = unit(0.5, "cm")))
ht_draw <- if (!is.null(row_mark)) row_mark + ht + rowRobustness else ht + rowRobustness
draw(ht_draw, heatmap_legend_side = "bottom")
dev.off()

pdf(file.path(dirOut, "TF_Robustness_Legend.pdf"), width = 2, height = 2)
draw(robustLegend)
dev.off()
saveRDS(k_clust, file.path(dirOut, "Kmeans_Jaccard_Edges_perTF.rds"))

# Aggregated TF Jaccard per network pair
for (fun_name in c("median", "mean")) {
  aggMat <- computeAggregatedPairJaccard(tfMatrix, fun = get(fun_name))
  plotAggregatedJaccard(aggMat,
                        file.path(dirOut, paste0("AggregatedTF_Jaccard_", fun_name, ".pdf")),
                        fig_width = 2.5, fig_height = 2.5)
}

# Hub TF overlap
for (metric in c("jaccard", "overlap")) {
  computeHubOverlapHeatmap(netDataList, tfCol = tfCol, networkNames = networkNames,
                           dirOut = dirOut, topN = 50, metric = metric,
                           fontsize = 9, fontsize_number = 7, fig_width = 2.5, fig_height = 2.5)
}

# Average/median targets per TF by robustness cluster
clust_df        <- data.frame(cluster = k_clust$cluster, TF = names(k_clust$cluster))
tfTargetCounts  <- tfTargetCounts %>% left_join(clust_df, by = "TF")
tfTargetCounts$avgTargets    <- rowMeans(tfTargetCounts[, networkNames], na.rm = TRUE)
tfTargetCounts$medianTargets <- apply(tfTargetCounts[, networkNames], 1, median, na.rm = TRUE)
write.table(tfTargetCounts,
            file.path(dirOut, "TF_avgTargetCounts_byRobustnessCluster.tsv"),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

for (yVar in c("avgTargets", "medianTargets")) {
  yLabel <- if (yVar == "avgTargets") "Average # of Targets" else "Median # of Targets"
  p <- ggplot(tfTargetCounts, aes(x = factor(cluster), y = .data[[yVar]],
                                   fill = factor(cluster))) +
    geom_boxplot(outlier.color = "red", outlier.size = 1, outlier.alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.1) +
    scale_fill_manual(values = annotCols) +
    labs(y = yLabel, x = "") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text.y  = element_text(size = 7, color = "black"),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.title   = element_text(size = 9, color = "black"),
          legend.position = "none",
          panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.1))
  ggsave(file.path(dirOut, paste0("TF_", yVar, "_Distribution.pdf")),
         plot = p, width = 2.5, height = 2, dpi = 600)
}

# ========================================================================
# PART 4: Degree distribution and stability plots
# Per network: targets per TF, TFs per target, stability of high/low degree genes
# Outputs saved alongside each network file
# ========================================================================

for (nm in networkNames) {
  data     <- netDataList[[nm]]
  basePath <- dirname(netFiles[[nm]])

  # Plot 1: Distribution of targets per TF
  dfTF <- data.frame(targetCount = tfTargetCounts[[nm]])
  p1 <- ggplot(dfTF, aes(x = targetCount)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "#0072B2", color = "black", alpha = 0.8) +
    labs(title = paste("Targets per TF —", nm), x = "# Targets", y = "# TFs") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 9),
          panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
          panel.grid.minor = element_blank())
  ggsave(file.path(basePath, paste0("TargetCountPerTF_", nm, ".pdf")),
         plot = p1, width = 3, height = 3, dpi = 600)

  # Plot 2: Distribution of TFs per target
  targetCounts   <- table(data[[targetCol]])
  dfTarget       <- data.frame(Target  = names(targetCounts),
                                tfCount = as.integer(targetCounts))
  dfTargetSorted <- dfTarget[order(dfTarget$tfCount), ]
  p2 <- ggplot(dfTarget, aes(x = tfCount)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "steelblue", color = "black", alpha = 0.8) +
    labs(title = paste("TFs per Target —", nm), x = "# TFs", y = "# Targets") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 9),
          panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
          panel.grid.minor = element_blank())
  ggsave(file.path(basePath, paste0("TFCountPerTarget_", nm, ".pdf")),
         plot = p2, width = 3, height = 3, dpi = 600)

  # Plot 3: Stability boxplot for top/bottom nSelect degree target genes
  lowTargs  <- head(dfTargetSorted$Target, nSelect)
  highTargs <- tail(dfTargetSorted$Target, nSelect)
  stabilityDF <- data[data[[targetCol]] %in% c(lowTargs, highTargs), c(targetCol, stabCol)]
  stabilityDF$Group <- ifelse(stabilityDF[[targetCol]] %in% lowTargs,
                               "Low TFs per Target", "High TFs per Target")
  stabilityDF <- merge(stabilityDF, dfTargetSorted, by.x = targetCol, by.y = "Target")
  stabilityDF$targetLabel <- factor(
    paste0(stabilityDF[[targetCol]], " (", stabilityDF$tfCount, ")"),
    levels = unique(paste0(stabilityDF[[targetCol]], " (", stabilityDF$tfCount, ")")[
      order(stabilityDF$tfCount)])
  )
  p3 <- ggplot(stabilityDF, aes(x = targetLabel, y = Stability, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    labs(title = paste("Per-Target Stability —", nm),
         x = paste0("Top/bottom ", nSelect, " targets by TF count"), y = "Stability") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text.x  = element_text(size = 7, angle = 90, vjust = 0.5),
          axis.text.y  = element_text(size = 7),
          axis.title   = element_text(size = 9),
          legend.position = "top",
          panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
          panel.grid.minor = element_blank())
  ggsave(file.path(basePath, paste0("Top", nSelect, "HighLow_inDegree_Stability_", nm, ".pdf")),
         plot = p3, width = 10, height = 5, dpi = 600)
}
