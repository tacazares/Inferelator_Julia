# =============================================================================
# DownSamplingObject.R
#
# Downsample a Seurat object to specified cell counts while preserving
# cell type proportions. Saves both raw counts and log-normalized counts
# as Arrow files for use with InferelatorJL.
#
# Two output files per downsampled size:
#   *_rawCounts.arrow      — raw integer counts  → use with pseudobulk + VST
#   *_logNorm.arrow        — log-normalized       → use with direct normalization
#
# Output format (both files):
#   - Rows    = genes
#   - Columns = cells
#   - First column named "Genes"
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
})

# =============================================================================
# USER INPUTS — edit this section
# =============================================================================

rdsPath     <- "/path/to/seurat_object.rds"
dirOut      <- "/path/to/output"
datasetName <- "dataset"              # used in output filenames
celltypeCol <- "CellType"             # metadata column with cell type labels
nList       <- c(30000, 10000, 1000)  # cell counts to downsample to
seed        <- 123

# =============================================================================

dir.create(dirOut, showWarnings = FALSE, recursive = TRUE)

obj      <- readRDS(rdsPath)
allCells <- colnames(obj)
cellTypes <- obj@meta.data[[celltypeCol]]

# Proportional sampling weights (preserves cell type composition)
nCelltype  <- table(cellTypes)
pctCelltype <- prop.table(nCelltype)
indexProbs <- pctCelltype[as.character(cellTypes)]

# Helper: format n as label (e.g. 10000 → "10K", 1000 → "1K")
.cellLabel <- function(n) paste0(sub("\\..*", "", as.character(n / 1000)), "K")

# Helper: convert matrix to Arrow-ready data.frame (Genes column first)
.toArrowDF <- function(mat) {
  df       <- as.data.frame(as.matrix(mat))
  df$Genes <- rownames(df)
  df[, c("Genes", setdiff(colnames(df), "Genes"))]
}

for (n in nList) {
  cellLabel <- .cellLabel(n)
  cat("\n--- Downsampling to", n, "cells (", cellLabel, ") ---\n")

  set.seed(seed)
  sampledIdx   <- sample(seq_along(cellTypes), size = n,
                         prob = indexProbs, replace = FALSE)
  sampledCells <- allCells[sampledIdx]
  objSubset    <- subset(obj, cells = sampledCells)

  # Log-normalize (creates layer = 'data')
  objSubset <- NormalizeData(objSubset, verbose = FALSE)

  # --- Raw counts (layer = 'counts') → for pseudobulk + VST pipeline ---
  rawCounts <- GetAssayData(objSubset, layer = "counts")
  write_feather(.toArrowDF(rawCounts),
                file.path(dirOut, paste0(datasetName, "_", cellLabel,
                                         "Cells_rawCounts.arrow")))
  cat("  Saved raw counts →", paste0(datasetName, "_", cellLabel, "Cells_rawCounts.arrow"), "\n")

  # --- Log-normalized counts (layer = 'data') → for direct normalization ---
  logNorm <- GetAssayData(objSubset, layer = "data")
  write_feather(.toArrowDF(logNorm),
                file.path(dirOut, paste0(datasetName, "_", cellLabel,
                                         "Cells_logNorm.arrow")))
  cat("  Saved log-norm  →", paste0(datasetName, "_", cellLabel, "Cells_logNorm.arrow"), "\n")
}

cat("\nDone. All outputs saved to:", dirOut, "\n")
