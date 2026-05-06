# =============================================================================
#  evaluateNetworks.R — Pairwise network comparison and structural analysis
#
#  Compares multiple GRN inference runs (e.g. bStARS vs bEBIC variants) across
#  four dimensions:
#
#  1. Summary statistics  — edge counts, TF/target coverage, prior support
#  2. Global similarity   — pairwise Jaccard across all edges and top-N% thresholds
#  3. TF-centric Jaccard  — per-TF regulon consistency; robustness clusters
#  4. Degree distribution — how edges are distributed across TFs and target genes
#
#  All outputs are organised into numbered subdirectories under dirOut so the
#  folder structure mirrors the analysis narrative.
#
# =============================================================================

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
  library(scales)
})

source("/path/to/evaluation/R/evaluateNetUtils.R")


# =============================================================================
# USER INPUTS — edit this section before running
# =============================================================================

# Root output directory — subdirectories are created automatically
dirOut <- "/path/to/output/networkComparison"

# Named list of networks to compare: label => path to combined edge file
# Use combined_max_sp.tsv (sparse/long format) for cross-network comparison
netFiles <- list(
  "NetworkA" = "/path/to/networkA/Combined/combined_max_sp.tsv",
  "NetworkB" = "/path/to/networkB/Combined/combined_max_sp.tsv"
)

# Prior network file (sparse format: TF, Target, Weight columns)
priorFile <- "/path/to/priors/prior_sp.tsv"

# Column names in edge files
tfCol     <- "TF"             # TF column name
targetCol <- "Gene"           # Target gene column name
rankCol   <- "signedQuantile" # Edge ranking score column
stabCol   <- "Stability"      # Edge confidence/stability column

# Prior column names
priorTfCol     <- "TF"
priorTargetCol <- "Target"
priorWgtCol    <- "Weight"

# Analysis parameters
nSelect      <- 10   # Top/bottom N target genes by TF count for stability boxplot
N            <- NULL # Top-N% edges for comparisons (NULL = use all edges)
tfList       <- NULL # Restrict TF Jaccard to a subset of TFs (NULL = all TFs)
k_center     <- 6    # Number of k-means clusters for TF Jaccard heatmap
file_k_clust <- NULL # Path to saved k-means RDS to reuse clustering (NULL = recompute)

# Marker TFs to annotate in the TF Jaccard heatmap (set to c() to skip)
lineageTFs <- c("TF1", "TF2", "TF3")


# =============================================================================
# OUTPUT SUBDIRECTORIES
# Numbered to enforce reading order and mirror the analysis narrative
# =============================================================================

dir01 <- file.path(dirOut, "01_summary")
dir02 <- file.path(dirOut, "02_global_similarity")
dir03 <- file.path(dirOut, "03_TF_jaccard")
dir04 <- file.path(dirOut, "04_degree_distribution")
for (d in c(dir01, dir02, dir03, dir04)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# Close any stale graphics devices left open by a previous partial run.
# A stale pdf() device with the same output path as a new call will corrupt
# the output file; clearing here prevents that.
while (dev.cur() > 1) dev.off()


# =============================================================================
# LOAD DATA
# =============================================================================

priorDt   <- read.table(priorFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
priorData <- melt(priorDt, id.vars = "X", value.name = "Weight")
colnames(priorData) <- c("Target", "TF", "Weight")
priorData  <- priorData[priorData[[priorWgtCol]] != 0, ]
priorPairs <- unique(paste(priorData[[priorTfCol]], priorData[[priorTargetCol]], sep = "_"))

netDataList  <- lapply(netFiles, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
networkNames <- names(netDataList)
numNetworks  <- length(netDataList)

# Shared theme for all ggplot figures
# Major gridlines: visible reference; minor: subtle subdivision
gridTheme <- theme(
  panel.grid.major   = element_line(color = "grey75",  linewidth = 0.3),
  panel.grid.minor   = element_line(color = "#e8e8e8", linewidth = 0.1),
  panel.grid.major.x = element_line(color = "grey75",  linewidth = 0.3),
  panel.grid.minor.x = element_line(color = "#e8e8e8", linewidth = 0.1),
  panel.border       = element_rect(color = "grey50",  linewidth = 0.3, fill = NA),
  axis.ticks         = element_line(linewidth = 0.15,  color = "grey60"),
  axis.ticks.length  = unit(2, "pt")
)


# =============================================================================
# PART 1: Summary statistics
#
# Story: Basic sanity check — how many TFs, targets, and edges does each
#        network have? What fraction of edges are supported by the prior?
#        Establishes that networks are comparable in scale before deeper analysis.
#
# Output: 01_summary/network_summary_statistics.tsv
# =============================================================================

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
write.table(summaryTable,
            file.path(dir01, "network_summary_statistics.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")


# =============================================================================
# PART 2: Global pairwise similarity
#
# Story: How similar are the networks to each other overall?
#        The global Jaccard heatmap gives the headline number.
#        The line plot shows whether similarity holds at all confidence
#        thresholds or only at the very top — a flat line means agreement is
#        consistent; a rising line means networks converge only at high confidence.
#
# Outputs:
#   02_global_similarity/GlobalJaccard_heatmap_allEdges.pdf
#   02_global_similarity/JaccardSimilarity_acrossTopNpercent.pdf
#   02_global_similarity/pairwise_edge_intersection_counts.tsv
#   02_global_similarity/edges_unique_to_each_network.tsv
#   02_global_similarity/jaccard_similarity_byTopNpercent.tsv
# =============================================================================

# Per-TF target counts (computed here, reused in Parts 3 and 4)
tfTargetCounts <- data.frame(
  TF = unique(unlist(lapply(netDataList, function(x) unique(x[[tfCol]])))))
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
fileSuffix <- if (!is.null(N)) paste0("top", N, "pct") else "allEdges"

# Edges exclusive to each network (not shared with any other)
uniqCountsDf <- data.frame(
  network        = names(edgeSets),
  numUniqueEdges = sapply(seq_along(edgeSets), function(k) {
    length(setdiff(edgeSets[[k]], unlist(edgeSets[-k])))
  }),
  stringsAsFactors = FALSE
)
write.table(uniqCountsDf,
            file.path(dir02, paste0("edges_unique_to_each_network_", fileSuffix, ".tsv")),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

# Global Jaccard heatmap — one number per network pair
pairs <- computeJaccardMatrix(edgeSets)
plotSimilarityHeatmap(pairs[["pairJac"]],
                      fileOut    = file.path(dir02, paste0("GlobalJaccard_heatmap_", fileSuffix, ".pdf")),
                      fig_width  = 4, fig_height = 4)

# Pairwise intersection counts (supporting table for heatmap)
pairDf <- melt(pairs[["pairInt"]], varnames = c("network1", "network2"),
               value.name = "numSharedEdges", na.rm = TRUE)
write.table(pairDf,
            file.path(dir02, paste0("pairwise_edge_intersection_counts_", fileSuffix, ".tsv")),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

# Jaccard across top-N% thresholds — does similarity hold at all confidence levels?
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
  mutate(Network1 = as.character(Network1),
         Network2 = as.character(Network2),
         pairID   = ifelse(Network1 < Network2,
                           paste(Network1, Network2, sep = "~"),
                           paste(Network2, Network1, sep = "~"))) %>%
  group_by(TopNpercent, pairID) %>% slice(1) %>% ungroup()

write.table(jaccardDf,
            file.path(dir02, "jaccard_similarity_byTopNpercent.tsv"),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

uniquePairs     <- sort(unique(jaccardDf$pairID))
comparisonColor <- setNames(
  colorRampPalette(brewer.pal(min(length(uniquePairs) + 1, 9), "Set1"))(length(uniquePairs)),
  uniquePairs
)

p_jac <- ggplot(jaccardDf, aes(x = TopNpercent, y = Jaccard, color = pairID, group = pairID)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = comparisonColor) +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x     = "Top N% edges (by confidence score)",
       y     = "Jaccard similarity",
       title = "Pairwise network similarity across confidence thresholds") +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text    = element_text(size = 7, color = "black"),
        axis.title   = element_text(size = 9, color = "black"),
        legend.title = element_blank(),
        legend.text  = element_text(size = 7)) +
  gridTheme

ggsave(file.path(dir02, "JaccardSimilarity_acrossTopNpercent.pdf"),
       plot = p_jac, width = 7, height = 4, dpi = 600)


# =============================================================================
# PART 3: TF-centric Jaccard
#
# Story: Global Jaccard treats all edges equally. Here we ask per TF:
#        how consistently does each method predict the same targets?
#        High-Jaccard TFs are robustly predicted regardless of method —
#        good candidates for biological follow-up.
#        Low-Jaccard TFs are method-sensitive — interpret with caution.
#
#        K-means clusters group TFs by their cross-network consistency profile.
#        Hub TF overlap asks whether methods agree on which TFs are most
#        connected — biologically the most consequential comparison.
#        Jaccard penalises size differences between hub sets; the overlap
#        coefficient does not, making both complementary.
#
# Outputs:
#   03_TF_jaccard/TFJaccard_heatmap_byComparison.pdf
#   03_TF_jaccard/TFJaccard_heatmap_byComparison_robustnessLegend.pdf
#   03_TF_jaccard/TFJaccard_aggregated_median.pdf
#   03_TF_jaccard/HubTF_overlap_Jaccard_top50.pdf
#   03_TF_jaccard/HubTF_overlap_coefficient_top50.pdf
#   03_TF_jaccard/TF_medianTargetCount_byRobustnessCluster.pdf
#   03_TF_jaccard/TF_targetCounts_byRobustnessCluster.tsv
#   03_TF_jaccard/kmeans_TFJaccard_clustering.rds
# =============================================================================

tfMatrix     <- computeTFJaccard(netDataList, tfList)
tfMatrix[is.na(tfMatrix)] <- 0
tfRobustness <- rowMeans(tfMatrix, na.rm = TRUE)

# K-means clustering (reuse saved result if provided to ensure reproducibility)
k_clust  <- if (is.null(file_k_clust)) {
  kmeans(tfMatrix, centers = k_center, nstart = 20, iter.max = 50)
} else {
  readRDS(file_k_clust)
}
split_by <- factor(k_clust$cluster, levels = seq_len(k_center))

# Color setup
heatCols   <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
robustCols <- colorRamp2(seq(min(tfRobustness), max(tfRobustness), length.out = 100),
                         colorRampPalette(brewer.pal(9, "Blues"))(100))
nPairs     <- ncol(tfMatrix)
compColor  <- setNames(
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
    at        = which(rownames(tfMatrix) %in% lineageTFs),
    labels    = rownames(tfMatrix)[rownames(tfMatrix) %in% lineageTFs],
    labels_gp = gpar(fontsize = 6, fontface = "plain"),
    padding   = unit(0.3, "mm"), side = "left"
  ))
} else { NULL }

robustLegend <- Legend(
  title     = "TF Robustness",
  col_fun   = robustCols,
  title_gp  = gpar(fontsize = 8, fontface = "plain"),
  labels_gp = gpar(fontsize = 6)
)

# Main TF Jaccard heatmap — rows = TFs, columns = pairwise comparisons
# tryCatch/finally ensures dev.off() is always called even if draw() errors,
# preventing stale pdf devices that would corrupt subsequent output files.
tryCatch({
  pdf(file.path(dir03, "TFJaccard_heatmap_byComparison.pdf"),
      width = 5, height = 3, compress = TRUE)
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
                  direction      = "horizontal", title = "Jaccard",
                  title_position = "topcenter",
                  title_gp       = gpar(fontsize = 8, fontface = "plain"),
                  labels_gp      = gpar(fontsize = 6),
                  legend_width   = unit(2.5, "cm"),
                  legend_height  = unit(0.5, "cm")))
  ht_draw <- if (!is.null(row_mark)) row_mark + ht + rowRobustness else ht + rowRobustness
  draw(ht_draw, heatmap_legend_side = "bottom")
}, finally = invisible(dev.off()))

# Robustness colour scale legend (saved separately for figure assembly)
tryCatch({
  pdf(file.path(dir03, "TFJaccard_heatmap_byComparison_robustnessLegend.pdf"),
      width = 2, height = 2)
  draw(robustLegend)
}, finally = invisible(dev.off()))

saveRDS(k_clust, file.path(dir03, "kmeans_TFJaccard_clustering.rds"))

# Aggregated TF Jaccard — median per network pair (one number per TF pair per comparison)
# Summarises the heatmap into a single comparable matrix
aggMat <- computeAggregatedPairJaccard(tfMatrix, fun = median)
plotAggregatedJaccard(aggMat,
                      file.path(dir03, "TFJaccard_aggregated_median.pdf"),
                      fig_width = 5, fig_height = 5)

# Hub TF overlap — do methods agree on the most connected TFs?
# Jaccard: penalises size differences between hub sets
# Overlap coefficient: size-insensitive — how much of the smaller set is captured?
# Both are complementary: use Jaccard if hub set sizes are similar,
# overlap coefficient if one method predicts systematically more hub TFs
for (metric in c("jaccard", "overlap")) {
  metricLabel <- if (metric == "jaccard") "Jaccard" else "overlap_coefficient"
  computeHubOverlapHeatmap(netDataList, tfCol = tfCol, networkNames = networkNames,
                           dirOut    = dir03,
                           topN      = 50,
                           metric    = metric,
                           fontsize  = 9, fontsize_number = 7,
                           fig_width = 5, fig_height = 5)
  # Note: computeHubOverlapHeatmap saves its own file; rename for clarity if needed
}

# Median targets per TF by robustness cluster
# Story: do high-robustness TFs tend to be hub TFs (many targets)?
clust_df       <- data.frame(cluster = k_clust$cluster, TF = names(k_clust$cluster))
tfTargetCounts <- tfTargetCounts %>% left_join(clust_df, by = "TF")
tfTargetCounts$medianTargets <- apply(tfTargetCounts[, networkNames], 1, median, na.rm = TRUE)

write.table(tfTargetCounts,
            file.path(dir03, "TF_targetCounts_byRobustnessCluster.tsv"),
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

p_clust <- ggplot(tfTargetCounts,
                  aes(x = factor(cluster), y = medianTargets, fill = factor(cluster))) +
  geom_boxplot(outlier.color = "red", outlier.size = 1, outlier.alpha = 0.7,
               linewidth = 0.3) +
  geom_jitter(width = 0.2, size = 0.1) +
  scale_fill_manual(values = annotCols) +
  labs(y = "Median # targets across networks",
       x = "Robustness cluster",
       title = "Target count by TF robustness cluster") +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text.y  = element_text(size = 7, color = "black"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title   = element_text(size = 9, color = "black"),
        legend.position = "none") +
  gridTheme

ggsave(file.path(dir03, "TF_medianTargetCount_byRobustnessCluster.pdf"),
       plot = p_clust, width = 2.5, height = 2, dpi = 600)


# =============================================================================
# PART 4: Degree distribution
#
# Story: How are edges distributed across TFs and target genes?
#
#   TF degree (targets per TF): shows whether a method concentrates edges on
#     a few dominant hub TFs or spreads them more evenly. On a log scale,
#     a curve shifted right in the ECDF tail means more high-degree hub TFs.
#     The ECDF (empirical CDF) is used instead of a boxplot because the
#     distribution is heavy-tailed — the biologically interesting question
#     lives in the tail, which a boxplot compresses into outlier dots.
#
#   Gene degree (TFs per gene): constrained by meanEdgesPerGene so distributions
#     are expected to be nearly identical across networks. Shown as a sanity
#     check that the edge cap is applied consistently.
#
#   Stability by gene degree: do target genes with many predicted regulators
#     also have higher-confidence individual edges, or are edges spread thin?
#
# Outputs:
#   04_degree_distribution/TFDegree_ECDF_logScale_allNetworks.pdf
#   04_degree_distribution/TargetCountPerTF_<network>.pdf       (one per network)
#   04_degree_distribution/TFsPerTarget_distribution_<network>.pdf
#   04_degree_distribution/Stability_highLow_geneDegree_<network>.pdf
# =============================================================================

for (nm in networkNames) {
  data <- netDataList[[nm]]

  # -- Plot 1: Distribution of targets per TF (log scale) --
  # Heavy-tailed: log scale reveals the shape; binwidth=0.1 on log10 = ~1.26x bins
  dfTF <- data.frame(targetCount = tfTargetCounts[[nm]])
  dfTF <- dfTF[dfTF$targetCount > 0, , drop = FALSE]   # exclude TFs with no edges in this network

  p1 <- ggplot(dfTF, aes(x = targetCount)) +
    geom_histogram(binwidth = 0.1, boundary = 0,
                   fill = "#0072B2", color = "#005a8e", linewidth = 0.2, alpha = 0.8) +
    scale_x_log10(labels       = scales::comma,
                  breaks       = c(1, 10, 100, 1000),
                  minor_breaks = c(2:9, seq(20, 90, 10), seq(200, 900, 100))) +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    annotation_logticks(sides     = "b",
                        linewidth = 0.2,
                        color     = "grey60",
                        short     = unit(1.5, "pt"),
                        mid       = unit(2.5, "pt"),
                        long      = unit(3.5, "pt")) +
    labs(title = paste("Targets per TF —", nm),
         x     = "# Targets (log scale)",
         y     = "# TFs") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text  = element_text(size = 7),
          axis.title = element_text(size = 9)) +
    gridTheme

  ggsave(file.path(dir04, paste0("TargetCountPerTF_", nm, ".pdf")),
         plot = p1, width = 3, height = 3, dpi = 600)

  # -- Plot 2: Distribution of TFs per target gene --
  # Expected to be similar across networks (controlled by meanEdgesPerGene).
  # binwidth=1 appropriate since range is ~1–20 integers.
  targetCounts   <- table(data[[targetCol]])
  dfTarget       <- data.frame(Target  = names(targetCounts),
                               tfCount = as.integer(targetCounts))
  dfTargetSorted <- dfTarget[order(dfTarget$tfCount), ]

  p2 <- ggplot(dfTarget, aes(x = tfCount)) +
    geom_histogram(binwidth = 1, boundary = 0.5,
                   fill = "#0072B2", color = "#005a8e", linewidth = 0.2, alpha = 0.8) +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    labs(title = paste("TFs per target gene —", nm),
         x     = "# TFs",
         y     = "# Target genes") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text  = element_text(size = 7),
          axis.title = element_text(size = 9)) +
    gridTheme

  ggsave(file.path(dir04, paste0("TFsPerTarget_distribution_", nm, ".pdf")),
         plot = p2, width = 3, height = 3, dpi = 600)

  # -- Plot 3: Edge stability for high vs low in-degree target genes --
  # Do genes with many TFs have higher-confidence individual edges?
  lowTargs    <- head(dfTargetSorted$Target, nSelect)
  highTargs   <- tail(dfTargetSorted$Target, nSelect)
  stabilityDF <- data[data[[targetCol]] %in% c(lowTargs, highTargs),
                      c(targetCol, stabCol)]
  stabilityDF$Group <- ifelse(stabilityDF[[targetCol]] %in% lowTargs,
                              "Low TFs per target", "High TFs per target")
  stabilityDF <- merge(stabilityDF, dfTargetSorted, by.x = targetCol, by.y = "Target")
  stabilityDF$targetLabel <- factor(
    paste0(stabilityDF[[targetCol]], " (", stabilityDF$tfCount, ")"),
    levels = unique(paste0(stabilityDF[[targetCol]], " (", stabilityDF$tfCount, ")")[
      order(stabilityDF$tfCount)])
  )

  p3 <- ggplot(stabilityDF, aes(x = targetLabel, y = Stability, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8, linewidth = 0.3) +
    labs(title = paste("Edge stability — high vs low in-degree targets —", nm),
         x     = paste0("Top/bottom ", nSelect, " target genes by # TFs"),
         y     = "Edge stability") +
    theme_bw(base_family = "Helvetica") +
    theme(axis.text.x     = element_text(size = 7, angle = 90, vjust = 0.5),
          axis.text.y     = element_text(size = 7),
          axis.title      = element_text(size = 9),
          legend.position = "top",
          legend.text     = element_text(size = 7),
          legend.title    = element_blank()) +
    gridTheme

  ggsave(file.path(dir04, paste0("Stability_highLow_geneDegree_", nm, ".pdf")),
         plot = p3, width = 10, height = 5, dpi = 600)
}

# -- Cross-network TF degree ECDF --
# Story: on a log scale, a curve shifted right in the tail means that network
#        produces more high-degree hub TFs. Overlapping curves indicate the
#        methods agree on how edges are distributed across TFs.
# ECDF is preferred over boxplot here because if a distribution is heavy-tailed —
# the hub TF differences live in the tail, which a boxplot compresses.
tfDegreeCombined <- do.call(rbind, lapply(networkNames, function(nm) {
  data.frame(network     = nm,
             targetCount = tfTargetCounts[[nm]])
}))
tfDegreeCombined <- tfDegreeCombined[tfDegreeCombined$targetCount > 0, ]

p_ecdf <- ggplot(tfDegreeCombined, aes(x = targetCount, color = network)) +
  stat_ecdf(linewidth = 0.7) +
  scale_x_log10(labels       = scales::comma,
                breaks       = c(1, 10, 100, 1000),
                minor_breaks = c(2:9, seq(20, 90, 10), seq(200, 900, 100))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent) +
  annotation_logticks(sides     = "b",
                      linewidth = 0.2,
                      color     = "grey60",
                      short     = unit(1.5, "pt"),
                      mid       = unit(2.5, "pt"),
                      long      = unit(3.5, "pt")) +
  labs(x     = "# Targets per TF (log scale)",
       y     = "Cumulative proportion of TFs",
       title = "TF degree distribution across networks") +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text       = element_text(size = 7),
        axis.title      = element_text(size = 9),
        legend.position = "bottom",
        legend.text     = element_text(size = 7),
        legend.title    = element_blank()) +
  gridTheme

ggsave(file.path(dir04, "TFDegree_ECDF_logScale_allNetworks.pdf"),
       plot = p_ecdf, width = 5, height = 4, dpi = 600)
