rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  # library(VennDiagram)
  library(ggVennDiagram)
  library(grid)
  library(reshape2)
})

source("/data/miraldiNB/Michael/Scripts/VennDiagram.R")

# dirOut <- "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/Bulk/5kbTSS/newPipe/spikeIgnored/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63"
dirOut <- "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/networkEval/"
dir.create(dirOut, recursive = T, showWarnings = F)

# Combined/combined_max.tsv
# TFmRNA/edges_subset.txt
# netFiles <- list(
#                 "TFA" =    "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/TFA/edges_subset.tsv",
#                 "TFmRNA" = "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/TFmRNA/edges_subset.tsv",
#                 "Combined" = "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/SC/SCENIC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5_eRegulon/Combined/combined_max.tsv"
#                 )

netFiles <- list(
                "TFA" =    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.tsv",
                "TFmRNA" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.tsv",
                "Combined" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/combined_max.tsv"
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
compareNets <- TRUE
# ========================================================================
# ------- Read Input files
# ======================================================================== 

# ----- Read Prior
priorData <- read.table(priorFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# priorData <- priorData[priorData[[priorWgtCol]] != 0, ]
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
if (compareNets){

  #--- Build edge sets
  edgeSets <- lapply(netDataList, function(df) paste(df[[tfCol]], df[[targetCol]], sep="~"))
  names(edgeSets) <- names(netDataList)
  networkNames <- names(netDataList)
  numNetworks <- length(netDataList)

  #--- Edges shared between all networks
  allShared <- Reduce(intersect, edgeSets)

  #1. --- Pairwise Jaccard, Intersections, and Heatmap

  if (numNetworks > 1) {
      pairInt  <- matrix(NA, numNetworks, numNetworks, dimnames=list(networkNames,networkNames))
      pairJac  <- matrix(NA, numNetworks, numNetworks, dimnames=list(networkNames,networkNames))

      for (i in 1:(numNetworks-1)) {
        for (j in (i+1):numNetworks) {
          shared <- intersect(edgeSets[[i]], edgeSets[[j]])
          un     <- union(edgeSets[[i]], edgeSets[[j]])
          pairInt[i,j] <- length(shared)
          pairJac[i,j] <- pairJac[j,i] <- if (length(un)) length(shared)/length(un) else NA
        }
      }
      diag(pairInt) <- vapply(edgeSets, length, integer(1))
      diag(pairJac) <- 1

      # heat-map of Jaccard ---------------------------------------- 
      heatCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(100)
      pdf(file.path(dirOut, "Jaccard.pdf"))
      pheatmap(pairJac, display_numbers = TRUE, main = "Pairwise Jaccard Index",  color = heatCols)
      dev.off()

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

      # 3. ------- Save Pairwise and uniqueDF counts
      pairDf <- melt(pairInt, varnames = c("network1", "network2"), value.name = "numIntersection", na.rm = TRUE)  
      write.table(pairDf, file.path(dirOut, "num_Pairwise_Intersection.tsv"), quote = F, col.names = NA, row.names=TRUE, sep = "\t")

  }

  # 2.------ Venn or UpSet Plot
  if (numNetworks >= 2) {
    vennInput <- edgeSets
    ggVennDiagram(vennInput, label_alpha = 0.6) + scale_fill_gradient(low="grey90", high = "red")
    ggsave(file.path(dirOut, "Venn2.pdf"))
    # venn.plot <- venn.diagram(
    #   vennInput,
    #   category.names = names(vennInput),
    #   filename = NULL,
    #   output = TRUE,
    #   main = "Shared TF~Gene Pairs"
    # )
    # pdf(file.path(dirOut, "Venn2.pdf"))
    # grid::grid.draw(venn.plot)
    # dev.off()

  } else if (numNetworks > 3) {
    # allEdges <- unique(unlist(edgeSets))
    # incMat <- sapply(edgeSets, function(edgeSet) as.integer(allEdges %in% edgeSet))
    # colnames(incMat) <- names(edgeSets)
    # rownames(incMat) <- allEdges
    upSetData <- edgeSets 

    plotUpSet <- function(setsData, Mode){
      mObj <- make_comb_mat(setsData, mode = Mode)
      ht <- UpSet(mObj, set_order = names(setsData),
                  top_annotation = upset_top_annotation(mObj, annotation_name_rot = 90, axis = FALSE,
                                          add_numbers = TRUE,  numbers_rot  = 90, 
                                          gp = gpar(col = comb_degree(mObj), fontsize = 6), height = unit(4, "cm"),
                                          axis_param = list(side = "left")),
                  right_annotation = upset_right_annotation(mObj,  axis_param = list(side = "bottom"), #labels = FALSE,labels_rot = 0
                                          gp = gpar(fill = "black", fontsize = 6),  width = unit(4, "cm"), show_annotation_name = FALSE,
                                          add_numbers = TRUE, axis = FALSE
                                                )
                )
      return(ht)
    }

    pdf(file.path(dirOut, "UpSet_distinct.pdf"), width = 8.5, height = 4)
    htDistinct <- plotUpSet(upSetData, Mode = "distinct")
    draw(htDistinct)
    dev.off()

    pdf(file.path(dirOut, "UpSet_Intersect.pdf"), width = 8.5, height = 4)
    htIntersect <- plotUpSet(upSetData, Mode = "intersect")
    draw(htIntersect)
    dev.off()
  }
}


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
  tfTargetCounts <- table(data[[tfCol]])
  dfTF <- data.frame(TF = names(tfTargetCounts), targetCount = as.integer(tfTargetCounts))

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

