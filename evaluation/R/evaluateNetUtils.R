suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(pheatmap)
  # library(VennDiagram)
  # library(ggVennDiagram)
  library(grid)
  library(reshape2)
  library(RColorBrewer)

})

# ==========================
# PART A: EDGE SET CONSTRUCTION
# ==========================

# -------- 1. Build Edge Sets -----------
getEdgeSets <- function(netDataList, tfCol="TF", targetCol="Gene", rankCol="signedQuantile", 
                        mode = c("topNpercent", "topN"), N = NULL) {
  #' getEdgeSets
  #'
  #' Extract edge sets from a list of networks, optionally selecting top-ranked edges.
  #'
  #' @param netDataList A named list of data.frames. Each data.frame represents a network with at least
  #'        the columns for TF, target gene, and ranking.
  #' @param tfCol Character. Name of the column containing transcription factors (default "TF").
  #' @param targetCol Character. Name of the column containing target genes (default "Gene").
  #' @param rankCol Character. Name of the column containing ranking scores (default "signedQuantile").
  #' @param mode Character. Method to select edges when N is specified. Choices are:
  #'        "topNpercent" (default) - select top N percent of edges based on rankCol,
  #'        "topN" - select top N edges.
  #'        Ignored if N is NULL.
  #' @param N Numeric. Number of edges (for "topN") or percent (for "topNpercent") to select.
  #'        If NULL (default), all edges are returned and mode is ignored.
  #'
  #' @return A named list of character vectors. Each element contains edges in the form "TF~Gene".
  #'
  mode <- match.arg(mode)
  
  edgeSets <- lapply(seq_along(netDataList), function(i) {
    df <- netDataList[[i]] %>% 
          filter(.data[[rankCol]] != 0) %>% 
          arrange(desc(.data[[rankCol]]))
    
    if(!rankCol %in% colnames(df)) stop(paste("Column", rankCol, "not found in network", names(netDataList)[i]))
    
    df <- df %>%
      mutate(edge = paste(.data[[tfCol]], .data[[targetCol]], sep="~"))
    # If N is NULL, return all edges
    if (is.null(N)) {
      return(df$edge)
    }
    
    if(mode == "topN") {
      df <- df %>% slice_head(n = N)
    } else { # topNpercent
      cutoff <- quantile(df[[rankCol]], 1 - N/100, na.rm = TRUE)
      df <- df %>% filter(abs(.data[[rankCol]]) >= cutoff) 
    }

    df$edge
  })
  
  names(edgeSets) <- names(netDataList)
  return(edgeSets)
}

# ======================================================
# PART B: Pairwise Overlap / Set Metrics & Correlation
# ======================================================

# - 1. Global Network Jaccard -----------
computeJaccardMatrix <- function(edgeSets) {
  networkNames <- names(edgeSets)
  numNetworks <- length(edgeSets)
  
  pairInt <- matrix(NA, nrow=numNetworks, ncol=numNetworks, dimnames=list(networkNames, networkNames))
  pairJac <- matrix(NA, nrow=numNetworks, ncol=numNetworks, dimnames=list(networkNames, networkNames))
  
  if (numNetworks > 1) {
    for (i in 1:(numNetworks-1)) {
      for (j in (i+1):numNetworks) {
        shared <- intersect(edgeSets[[i]], edgeSets[[j]])
        un     <- union(edgeSets[[i]], edgeSets[[j]])
        pairInt[i,j] <- length(shared)
        pairJac[i,j] <- pairJac[j,i] <- if (length(un)) length(shared)/length(un) else NA
      }
    }
  }
  
  diag(pairInt) <- vapply(edgeSets, length, integer(1))
  diag(pairJac) <- 1
  
  return(list(pairJac = pairJac, pairInt = pairInt))
}


# ==============================================================================
# - 2. Compute Spearman Correlation between pairs of ranks or weights
# ==============================================================================
#'
#' @param netDataList Named list of data.frames. Each data.frame must contain TF, target gene, and score columns.
#' @param edgeSets Optional named list of character vectors (from getEdgeSets). If provided, only these edges are used.
#' @param tfCol Character. Name of TF column (default "TF").
#' @param targetCol Character. Name of target gene column (default "Gene").
#' @param rankCol Character. Name of score column (default "signedQuantile").
#'
#' @return A tibble with columns: Net1, Net2, Spearman correlation.
#'
#' @examples
#' # Spearman on all edges
#' computeSpearman(netDataList)
#'
#' # Spearman on top edges only
#' edges_top10pct <- getEdgeSets(netDataList, N=10)
#' computeSpearman(netDataList, edgeSets = edges_top10pct)
computeSpearman <- function(netDataList, 
                            edgeSets = NULL,
                            ref = NULL,
                            tfCol="TF", 
                            targetCol="Gene", 
                            rankCol="signedQuantile") {
  
  prepForCor <- function(df, edges = NULL) {
    df <- df %>%
      filter(.data[[rankCol]] != 0) %>% 
      mutate(edge = paste(.data[[tfCol]], .data[[targetCol]], sep="~")) %>%
      dplyr::select(edge, score = .data[[rankCol]])
    
    if(!is.null(edges)) {
      df <- df %>% filter(edge %in% edges)
    }
    
    return(df)
  }
  
  # Prepare each network
  nets <- lapply(seq_along(netDataList), function(i) {
    edges_to_use <- if(!is.null(edgeSets)) edgeSets[[names(netDataList)[i]]] else NULL
    prepForCor(netDataList[[i]], edges = edges_to_use)
  })

  netNames <- names(netDataList)
  names(nets) <- netNames

  # ---  Compute pairwise Spearman correlations
  # Build comparison pairs
  if(!is.null(ref)){
    if(!ref %in% netNames){
      stop("Reference not found in edgeSets")
    }
    combs <- lapply(setdiff(netNames, ref), function(x)
      c(ref, x))
  } else {
    combs <- combn(netNames, 2, simplify = FALSE)
  }
  
  map_dfr(combs, function(p) {
    a <- nets[[p[1]]]
    b <- nets[[p[2]]]
    
    joined <- inner_join(a, b, by = "edge", suffix = c(".1", ".2"))
    rho <- if(nrow(joined) > 1) {
      cor(joined$score.1, joined$score.2, method = "spearman", use = "complete.obs")
    } else {NA}

    tibble(
      Net1 = p[1],
      Net2 = p[2],
      Spearman = rho
    )
  })
}


# ==== 3.Compute Edge  overlaps ====
computeOverlapStats <- function(edgeSets, ref = NULL){
  nets <- names(edgeSets)
  # Build comparison pairs
  if(!is.null(ref)){
    if(!ref %in% nets){
      stop("Reference not found in edgeSets")
    }
    combs <- lapply(setdiff(nets, ref), function(x)
      c(ref, x))
  } else {
    combs <- combn(nets, 2, simplify = FALSE)
  }
  # combs <- combn(names(edgeSets), 2, simplify=FALSE)
  
  map_dfr(combs, function(p){
    A <- edgeSets[[p[1]]]
    B <- edgeSets[[p[2]]]
    
    nA <- length(A)
    nB <- length(B)
    
    shared <- length(intersect(A,B))
    un  <- length(union(A,B))
    
    # Dice Coefficient (or Sorensen-Dice Index)
    dice <- ifelse(
      nA + nB > 0,
      2 * shared / (nA + nB),
      NA
    )
    
    # What fraction of the smaller set is shared with the larger set?
    # Szymkiewicz–Simpson coefficient
    fracOverlap =
      ifelse(min(nA,nB)>0,
             shared / min(nA,nB),
             NA)
    tibble(
      Net1 = p[1],
      Net2 = p[2],
      pairID = paste0(p[1], "~", p[2]),
      size1 = nA,
      size2 = nB,
      Intersection = shared,
      frac1in2 = ifelse(nA > 0, shared / nA, NA),
      frac2in1 = ifelse(nB > 0, shared / nB, NA),
      Jaccard = ifelse(un> 0, shared/un, NA),
      Dice = dice,
      FractionOverlap = fracOverlap
      
    )
  })
}

# ==== 3.Compute TF overlaps ====
computeTFoverlap <- function(netDataList, ref = NULL, tfCol = "TF") {
  # Build TF sets
  tfSets <- lapply(netDataList, function(df) unique(df[[tfCol]]))
  
  nets <- names(tfSets)
  
  # Build comparison pairs
  if(!is.null(ref)){
    if(!ref %in% nets){
      stop("Reference not found in netDataList")
    }
    combs <- lapply(setdiff(nets, ref), function(x) c(ref, x))
  } else {
    combs <- combn(nets, 2, simplify = FALSE)
  }
  
  # Compute pairwise overlaps
  map_dfr(combs, function(p){
    A <- tfSets[[p[1]]]
    B <- tfSets[[p[2]]]
    nA <- length(A)
    nB <- length(B)
    ov <- intersect(A,B)
    shared <- length(ov)
    un <- length(union(A,B))
    
    # dice <- ifelse(nA + nB > 0, 2 * shared / (nA + nB), NA)
    # fracOverlap <- ifelse(min(nA,nB) > 0, shared / min(nA,nB), NA)
    
    tibble(
      Net1 = p[1],
      Net2 = p[2],
      pairID = paste0(p[1], "~", p[2]),
      size1 = nA,
      size2 = nB,
      Overlaps = paste(ov, collapse = ","),
      nIntersection = shared,
      frac1in2 = ifelse(nA > 0, shared / nA, NA),
      frac2in1 = ifelse(nB > 0, shared / nB, NA)
      # Jaccard = ifelse(un > 0, shared/un, NA),
      # Dice = dice,
      # FractionOverlap = fracOverlap
    )
  })
}

# ====================================================
# PART C: TF-Level / Aggregated Metrics
# ====================================================

# ===== 1. TF-centric Jaccard ====
#' For each TF, how consistent are its predicted targets across networks?
#' Biological networks are TF-driven, so checking target consistency per TF
#' Jaccard penalizes sets with different sizes. It measures proportion of shared elements out of all unique elements
#' 
computeTFJaccard <- function(netDataList, tfList=NULL){
  networkNames <- names(netDataList)
  numNetworks <- length(netDataList)
  
  allTFs <- unique(unlist(lapply(netDataList, function(df) df$TF)))
  tfs <- if(!is.null(tfList)) intersect(allTFs, tfList) else allTFs
  
  pairNames <- combn(networkNames,2,function(x) paste(x, collapse="~"))
  tfMatrix <- matrix(NA, nrow=length(tfs), ncol=length(pairNames), dimnames=list(tfs, pairNames))
  
  for(tf in tfs){
    tfTargets <- lapply(netDataList, function(df) df$Gene[df$TF == tf])
    for(k in seq_along(pairNames)){
      pair <- combn(seq_len(numNetworks),2)[,k]
      shared <- intersect(tfTargets[[pair[1]]], tfTargets[[pair[2]]])
      un <- union(tfTargets[[pair[1]]], tfTargets[[pair[2]]])
      tfMatrix[tf,k] <- if(length(un)) length(shared) /length(un) else NA
    }
  }
  return(tfMatrix)
}

# plotTFHeatmap <- function(tfMatrix, fileOut, fontsize_number=7, fontsize=9, fig_width = 6, fig_height = 6){
#   heatCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(100)
#   pdf(fileOut, width = fig_width, height = fig_height)
#   pheatmap(tfMatrix,
#           display_numbers=TRUE,
#           number_color="black",
#           number_format="%.2f",
#           fontsize_number=fontsize_number,
#           fontsize=fontsize,
#           cluster_rows=FALSE,
#           cluster_cols=FALSE,
#           show_rownames = TRUE,
#           show_colnames = TRUE,
#           legend = FALSE ,   # removes grid lines
#           color=heatCols,
#           main=""
#             )
#   dev.off()
# }

# ==== 2. Aggregated TF Jaccard per Network Pair =====
computeAggregatedPairJaccard <- function(tfMatrix, fun=median){
  # tfMatrix: rows = TFs, cols = network pairs
  aggVec <- apply(tfMatrix, 2, fun, na.rm=TRUE)
  
  # Convert to symmetric matrix
  pairNames <- colnames(tfMatrix)
  nets <- unique(unlist(strsplit(pairNames, "~")))
  numNetworks <- length(nets)
  aggMat <- matrix(NA, nrow=numNetworks, ncol=numNetworks, dimnames=list(nets,nets))
  
  for(k in seq_along(pairNames)){
    pair <- strsplit(pairNames[k], "~")[[1]]
    aggMat[pair[1], pair[2]] <- aggVec[k]
    aggMat[pair[2], pair[1]] <- aggVec[k]
  }
  diag(aggMat) <- 1
  return(aggMat)
}


#' ==== 3./ Compute and plot hub TF set overlap across networks ====
#'
#' This function identifies the top N hub transcription factors (TFs) in each
#' network, computes pairwise similarity between hub sets across networks using
#' either the Jaccard index or the overlap coefficient, and plots a heatmap of
#' the similarity matrix.
#'
#' @param netDataList List of data frames, one per network, each containing edges.
#' @param tfCol Column name (string) indicating which column contains TF identifiers.
#' @param networkNames Character vector of network names, matching the order of netDataList.
#' @param dirOut Output directory for saving the heatmap PDF.
#' @param topN Integer. Number of top hubs (by edge count) to include per network. Default = 50.
#' @param metric Similarity metric. One of "jaccard" or "overlap". Default = "jaccard".
#' @param heatCols Color palette for the heatmap. Default = 100-color YlOrRd.
#' @param fontsize Base font size for heatmap text. Default = 9.
#' @param fontsize_number Font size for numbers displayed in cells. Default = 7.
#' #' @param fig_width Width of the PDF figure in inches. Default = 6.
#' @param fig_height Height of the PDF figure in inches. Default = 6.

#'
#' @return Invisibly returns the similarity matrix used to generate the heatmap.
#'
#' @examples
#' hubSim <- computeHubOverlapHeatmap(
#'   netDataList = netDataList,
#'   tfCol = "TF",
#'   networkNames = c("Net1", "Net2", "Net3"),
#'   dirOut = "results/",
#'   topN = 50,
#'   metric = "jaccard",
#'   fig_width = 8,
#'   fig_height = 8
#' )
#'
#' @export
computeHubOverlapHeatmap <- function(netDataList, tfCol, networkNames, 
                                     dirOut, topN = 50, 
                                     metric = c("jaccard", "overlap"),
                                     heatCols = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
                                     fontsize = 9, fontsize_number = 7, fig_width = 6, fig_height = 6) {
  # --- Argument check
  metric <- match.arg(metric)
  numNetworks <- length(netDataList)
  
  # --- 1. Count edges per TF and rank hubs
  hubRankList <- lapply(netDataList, function(df) {
    tfCounts <- table(df[[tfCol]])
    sort(tfCounts, decreasing = TRUE)
  })
  
  # --- 2. Take top N hubs
  topHubs <- lapply(hubRankList, function(x) head(names(x), topN))
  
  # --- 3. Initialize similarity matrix
  hubSim <- matrix(NA, numNetworks, numNetworks,
                   dimnames = list(networkNames, networkNames))
  
  # --- 4. Compute similarity
  for (i in 1:(numNetworks-1)) {
    for (j in (i+1):numNetworks) {
      inter <- intersect(topHubs[[i]], topHubs[[j]])
      
      if(metric == "jaccard") {
        union <- union(topHubs[[i]], topHubs[[j]])
        hubSim[i,j] <- hubSim[j,i] <- length(inter) / length(union)
      } else if(metric == "overlap") {
        minSize <- min(length(topHubs[[i]]), length(topHubs[[j]]))
        hubSim[i,j] <- hubSim[j,i] <- length(inter) / minSize
      }
    }
  }
  diag(hubSim) <- 1
  
  # --- 5. Plot heatmap
  pdf(file.path(dirOut, paste0("HubTF_Overlap_", metric, "_top", topN, ".pdf")), width = fig_width, height = fig_height)
  pheatmap(hubSim,
           display_numbers = TRUE,
           number_color = "black",
           number_format = "%.2f",
           fontsize_number = fontsize_number,
           fontsize = fontsize,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           legend = FALSE,
           color = heatCols,
           main = "")
  dev.off()
  
  invisible(hubSim)
}

# ==================================================
# - Plotting Ulis
# ==================================================


# ==== Edge Sharing ====
plotEdgeSharing <- function(edgeSets, fileSuffix, outDir=".") {
  library(ggVennDiagram)
  library(ComplexUpset)
  library(ComplexHeatmap)
  
  numNetworks <- length(edgeSets)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  
  if(numNetworks <= 3){
    # --- Venn diagram
    pdf(file.path(outDir, paste0("Venn2_", fileSuffix, ".pdf")), width=6, height=6)
    ggVennDiagram(edgeSets, label_alpha = 0.6) + 
      scale_fill_gradient(low="grey90", high="red") +
      theme(legend.position="none")
    dev.off()
    
  } else {
    # --- UpSet plots for >2 networks
    # Convert edgeSets to a presence/absence matrix
    allEdges <- unique(unlist(lapply(edgeSets, function(df) paste(df$TF, df$Gene, sep="~"))))
    incMat <- sapply(edgeSets, function(df) as.integer(allEdges %in% paste(df$TF, df$Gene, sep="~")))
    colnames(incMat) <- names(edgeSets)
    rownames(incMat) <- allEdges
    
    # Intersect mode
    mObjIntersect <- make_comb_mat(incMat, mode="intersect")
    pdf(file.path(outDir, "UpSet_Intersect.pdf"), width=8.5, height=4)
    draw(UpSet(mObjIntersect,
               set_order = names(edgeSets),
               top_annotation = upset_top_annotation(mObjIntersect, annotation_name_rot=90, axis=FALSE,
                                                     add_numbers=TRUE, numbers_rot=90,
                                                     gp=gpar(col=comb_degree(mObjIntersect), fontsize=6), height=unit(4,"cm")),
               right_annotation = upset_right_annotation(mObjIntersect, gp=gpar(fill="black", fontsize=6),
                                                         width=unit(4,"cm"), show_annotation_name=FALSE, add_numbers=TRUE)))
    dev.off()
    
    # Distinct mode
    mObjDistinct <- make_comb_mat(incMat, mode="distinct")
    pdf(file.path(outDir, paste0("UpSet_distinct_", fileSuffix, ".pdf")), width=8.5, height=4)
    draw(UpSet(mObjDistinct,
               set_order = names(edgeSets),
               top_annotation = upset_top_annotation(mObjDistinct, annotation_name_rot=90, axis=FALSE,
                                                     add_numbers=TRUE, numbers_rot=90,
                                                     gp=gpar(col=comb_degree(mObjDistinct), fontsize=6), height=unit(4,"cm")),
               right_annotation = upset_right_annotation(mObjDistinct, gp=gpar(fill="black", fontsize=6),
                                                         width=unit(4,"cm"), show_annotation_name=FALSE, add_numbers=TRUE)))
    dev.off()
  }
  
  message("Edge sharing plots generated.")
}


plotSimilarityHeatmap <- function(simMat, fileOut="sim_heatmap.pdf", 
                               fontsize_number=7, fontsize=9, fig_width=6, fig_height=6) {
  simMat[lower.tri(simMat)] <- t(simMat)[lower.tri(simMat)]
  heatCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(100)
  pdf(fileOut, height = fig_height, width = fig_width)
  pheatmap::pheatmap(simMat,
                     display_numbers=TRUE,
                     number_color="black",
                     number_format="%.2f",
                     fontsize_number=fontsize_number,
                     fontsize=fontsize,
                     cluster_rows=FALSE,
                     cluster_cols=FALSE,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     legend = FALSE , 
                     color=heatCols,
                     main=""
  )
  dev.off()
}

plotAggregatedJaccard <- function(aggMat, fileOut, fontsize_number=7, fontsize=9, fig_width = 6, fig_height = 6){
  heatCols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(100)
  pdf(fileOut, width = fig_width, height = fig_height)
  pheatmap(aggMat,
           display_numbers=TRUE,
           number_color="black",
           number_format="%.2f",
           fontsize_number=fontsize_number,
           fontsize=fontsize,
           cluster_rows=FALSE,
           cluster_cols=FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           legend = FALSE ,   # removes grid lines
           color=heatCols,
           main="Aggregated TF Jaccard per Network Pair")
  dev.off()
}

plotMetricTrend <- function(df, xCol, yCol, groupCol = "pairID", colorMap = NULL,
                            xLabel = NULL, yLabel = NULL, yLimits = NULL, xBreaks = NULL,
                            outFile = NULL, width = 6, height = 3, dpi = 600){
  p <- ggplot(df,
          aes_string(x = xCol, y = yCol, color = groupCol, group = groupCol)
      ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    theme_bw(base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.title  = element_text(size = 9, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.box.just = "left",
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.1)
    )

  # Optional scales
  if(!is.null(colorMap)){
    p <- p + scale_color_manual(values = colorMap)
  }

  if(!is.null(xBreaks)){
    p <- p + scale_x_continuous(breaks = xBreaks)
  }

  if(!is.null(yLimits)){
      p <- p + scale_y_continuous(limits = yLimits)
  } else{
      dataMin <- min(df[[yCol]], na.rm = TRUE)
      dataMax <- max(df[[yCol]], na.rm = TRUE)
      yMin <- if (dataMin >= 0) 0 else dataMin
      yMax <- max(dataMax, 1)
      yLimits <- c(yMin, yMax)
      p <- p + scale_y_continuous(limits = yLimits)
  }
  # Labels
  p <- p + labs( x = xLabel, y = yLabel )
  # Save if requested
  if(!is.null(outFile)) ggsave(outFile, plot = p, width = width, height = height, dpi = dpi)

  return(p)
}

# ==========================================
#  High-Level Pipeline / API (Top Level)
# ==========================================
computeTopNMetrics <- function(netDataList, 
                               ref = NULL,
                               topNs = NULL, 
                               tfCol="TF", 
                               targetCol="Target", 
                               rankCol="Weight",
                               mode = c("topNpercent", "topN")) {
  #' Compute Top-N overlap and correlation metrics for multiple networks
  #' 
  #' @param netDataList Named list of network data.frames
  #' @param topNs Numeric vector of percentages (Top N%) or number to evaluate
  #' @param tfCol Character, column name for transcription factor
  #' @param targetCol Character, column name for target gene
  #' @param rankCol Character, column name for ranking score
  #' 
  #' @return data.frame with columns: pairID, Jaccard, FractionOverlap, Spearman, topNs
  mode = match.arg(mode)
  if (is.null(topNs)) topNs <- "all"
  topNResults <- lapply(topNs, function(k) {
    cat("Top-N:", k, "\n")
    
    if (k == "all"){
      edges <- getEdgeSets(netDataList, tfCol = tfCol, 
                          targetCol = targetCol, rankCol = rankCol)
    } else{
    # Extract Top-N% edges
      edges <- getEdgeSets(netDataList, tfCol = tfCol, 
                           targetCol = targetCol, rankCol = rankCol,
                           mode = mode, N = k)
    }
    # Overlap statistics
    stats <- computeOverlapStats(edges, ref)
    stats$Npct <- if(k == "all") NA else k
    
    # Spearman correlation
    corRes <- computeSpearman(netDataList, edgeSets = edges, ref = ref, tfCol = tfCol, targetCol = targetCol, rankCol = rankCol)
    # Merge Spearman with stats
    stats <- stats %>%
      left_join(corRes %>%
                  mutate(pairID = paste(Net1, Net2, sep="~")) %>%
                  dplyr::select(pairID, Spearman),
                by = "pairID")
    
    stats
  })
  
  bind_rows(topNResults)
}

