# Network Evaluation Scripts (R)

Standalone R scripts for evaluating and comparing GRN outputs from InferelatorJL.
These scripts are not part of the Julia package — they are post-hoc diagnostic tools
that operate on the edge files produced by the pipeline.

---

## Files

| File | Description |
|---|---|
| `evaluateNetUtils.R` | Reusable utility functions for network comparison (source this in your scripts) |
| `evaluateNetworks.R` | End-to-end usage example — edit paths and run |
| `histogramConfidences.R` | Standalone script to plot edge confidence/stability distributions |

---

## Requirements

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "pheatmap",
                   "reshape2", "RColorBrewer", "purrr"))

# Bioconductor
BiocManager::install(c("ComplexHeatmap"))

# Optional (for Venn/UpSet plots)
install.packages("ggVennDiagram")
BiocManager::install("ComplexUpset")
```

---

## Input format

All functions expect network edge files in **long format** (tab-delimited `.tsv`),
as produced by InferelatorJL's `buildNetwork` step:

| Column | Description |
|---|---|
| `TF` | Transcription factor name |
| `Gene` | Target gene name |
| `signedQuantile` | Signed ranking score (used for edge selection) |
| `Stability` | Bootstrap stability score (0–max_subsamples) |
| `inPrior` | Whether the edge is supported by the prior (0/1) |

The prior file should be in **sparse format** with columns `TF`, `Target`, `Weight`.

---

## Quickstart

1. Open `evaluateNetworks.R`
2. Edit the **USER INPUTS** section at the top:
   - Set `dirOut` — where outputs will be saved
   - Set `netFiles` — named list of your network edge files
   - Set `priorFile` — path to the prior network
   - Adjust column names if needed
3. Run the script

---

## What `evaluateNetworks.R` produces

### Part 1 — Summary statistics
| Output | Description |
|---|---|
| `summaryStatistics.tsv` | Per-network counts: unique TFs, targets, total interactions, % prior-supported |

### Part 2 — Pairwise network comparison
| Output | Description |
|---|---|
| `Global_Jaccard_*.pdf` | Heatmap of pairwise Jaccard similarity (all edges) |
| `num_Pairwise_Intersection_*.tsv` | Number of shared edges between each pair |
| `num_UniqueEdges_*.tsv` | Edges unique to each network |
| `Jaccard_Sim.tsv` | Jaccard across top 10%–100% confidence thresholds |
| `linePlot_JaccardSim.pdf` | Line plot of Jaccard vs confidence threshold |

### Part 3 — TF-centric Jaccard
Measures how consistently each TF's predicted targets are recovered across networks.
A high Jaccard means the TF's regulon is robust; a low Jaccard means it is network-specific.

| Output | Description |
|---|---|
| `Jaccard_Edges_perTF.pdf` | Heatmap of per-TF Jaccard with k-means clustering and robustness bar |
| `TF_Robustness_Legend.pdf` | Standalone legend for the robustness color bar |
| `AggregatedTF_Jaccard_median.pdf` | Median per-TF Jaccard aggregated per network pair |
| `AggregatedTF_Jaccard_mean.pdf` | Mean per-TF Jaccard aggregated per network pair |
| `HubTF_Overlap_jaccard_top50.pdf` | Jaccard similarity of top-50 hub TFs across networks |
| `HubTF_Overlap_overlap_top50.pdf` | Overlap coefficient of top-50 hub TFs across networks |
| `Kmeans_Jaccard_Edges_perTF.rds` | Saved k-means clustering result for reuse |

### Part 4 — Degree distributions (saved alongside each network file)
| Output | Description |
|---|---|
| `TargetCountPerTF_<network>.pdf` | Distribution of number of targets per TF |
| `TFCountPerTarget_<network>.pdf` | Distribution of number of TFs per target gene |
| `Top<N>HighLow_inDegree_Stability_<network>.pdf` | Stability boxplot for high/low degree target genes |

---

## Key functions in `evaluateNetUtils.R`

### Edge set construction
```r
# Extract edges from a list of networks; optionally select top N% or top N edges
edgeSets <- getEdgeSets(netDataList, tfCol="TF", targetCol="Gene",
                         rankCol="signedQuantile", mode="topNpercent", N=NULL)
```

### Pairwise metrics
```r
# Pairwise Jaccard similarity matrix
computeJaccardMatrix(edgeSets)

# Pairwise overlap statistics (Jaccard, Dice, FractionOverlap, Intersection)
computeOverlapStats(edgeSets, ref=NULL)

# Pairwise Spearman correlation of edge rankings
computeSpearman(netDataList, edgeSets=NULL, ref=NULL)

# High-level: compute all metrics across multiple top-N% thresholds
computeTopNMetrics(netDataList, topNs=seq(10,100,10))
```

### TF-centric metrics
```r
# Per-TF Jaccard: how consistent is each TF's regulon across networks?
tfMatrix <- computeTFJaccard(netDataList, tfList=NULL)

# Aggregate per-TF Jaccard to a single value per network pair
computeAggregatedPairJaccard(tfMatrix, fun=median)

# TF overlap statistics (similar to edge overlap but for TF sets)
computeTFoverlap(netDataList, ref=NULL)

# Hub TF overlap heatmap (top N most connected TFs)
computeHubOverlapHeatmap(netDataList, tfCol, networkNames, dirOut, topN=50)
```

### Plotting
```r
# Similarity heatmap (works for any square similarity matrix)
plotSimilarityHeatmap(simMat, fileOut, fig_width=6, fig_height=6)

# Aggregated TF Jaccard heatmap
plotAggregatedJaccard(aggMat, fileOut, fig_width=6, fig_height=6)

# Line plot of a metric across top-N thresholds
plotMetricTrend(df, xCol, yCol, groupCol="pairID", outFile=NULL)

# Venn (<=3 networks) or UpSet (>3 networks) edge sharing plot
plotEdgeSharing(edgeSets, fileSuffix, outDir=".")
```
