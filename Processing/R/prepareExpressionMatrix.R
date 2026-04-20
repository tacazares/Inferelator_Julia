# =============================================================================
# prepareExpressionMatrix.R
#
# Converts gene expression data from common formats into Apache Arrow (.arrow)
# for use as input to InferelatorJL.
#
# InferelatorJL expects an Arrow file where:
#   - Rows are genes
#   - Columns are samples/cells
#   - First column is named "Genes" and contains gene names
#   - Values are normalized expression (e.g. log-normalized counts)
#
# Supported input formats:
#   1. Seurat object (.rds)
#   2. SingleCellExperiment object (.rds)
#   3. Dense matrix TSV/CSV (genes x samples)
#   4. 10x Genomics sparse matrix (matrix.mtx + barcodes.tsv + features.tsv)
#   5. HDF5 (.h5) — 10x Cell Ranger output
#   6. AnnData (.h5ad) — scanpy/Python output
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow)
})

# =============================================================================
# USER INPUTS — edit this section
# =============================================================================

# Choose input format:
# "seurat", "sce", "tsv", "10x", "h5", "h5ad"
inputFormat <- "seurat"

# Input file or folder path
inputPath <- "/path/to/input"

# Output Arrow file path
outputPath <- "/path/to/output/expression_logNorm.arrow"

# --- Seurat / SCE options ---
# Layer/assay to extract (log-normalized counts)
# Seurat v5: typically "data" (log-normalized) or "counts" (raw)
# SCE: "logcounts" (log-normalized) or "counts" (raw)
assayLayer <- "data"       # for Seurat
assayName  <- "RNA"        # Seurat assay name (usually "RNA")
sceAssay   <- "logcounts"  # for SingleCellExperiment

# --- TSV/CSV options ---
sep       <- "\t"    # delimiter: "\t" for TSV, "," for CSV
geneCol   <- 1       # column index or name containing gene names (set to NULL if genes are rownames)
hasHeader <- TRUE    # whether the file has a header row

# --- 10x options ---
# inputPath should be the folder containing matrix.mtx, barcodes.tsv, features.tsv
# Set normalize10x = TRUE to log-normalize raw counts (recommended)
normalize10x <- TRUE

# --- HDF5 options ---
# inputPath should be the .h5 file from Cell Ranger
normalizeH5 <- TRUE

# --- AnnData options ---
# inputPath should be the .h5ad file
# layer: NULL = use X (default), or name of a layer e.g. "log_norm"
h5adLayer <- NULL

# =============================================================================
# HELPER: format matrix as data.frame with Genes as first column
# =============================================================================

.toArrowDF <- function(mat) {
  df       <- as.data.frame(as.matrix(mat))
  df$Genes <- rownames(df)
  df       <- df[, c("Genes", setdiff(colnames(df), "Genes"))]
  df
}

# =============================================================================
# CONVERT
# =============================================================================

df <- switch(inputFormat,

  # ---------------------------------------------------------------------------
  # 1. Seurat object
  # ---------------------------------------------------------------------------
  "seurat" = {
    suppressPackageStartupMessages(library(Seurat))
    obj        <- readRDS(inputPath)
    norm_counts <- GetAssayData(obj, assay = assayName, layer = assayLayer)
    .toArrowDF(norm_counts)
  },

  # ---------------------------------------------------------------------------
  # 2. SingleCellExperiment object
  # ---------------------------------------------------------------------------
  "sce" = {
    suppressPackageStartupMessages(library(SingleCellExperiment))
    obj        <- readRDS(inputPath)
    norm_counts <- assay(obj, sceAssay)
    .toArrowDF(norm_counts)
  },

  # ---------------------------------------------------------------------------
  # 3. Dense TSV / CSV (genes x samples)
  # ---------------------------------------------------------------------------
  "tsv" = {
    mat <- read.table(inputPath, header = hasHeader, sep = sep,
                      stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(geneCol)) {
      genes <- mat[[geneCol]]
      mat   <- mat[, -geneCol, drop = FALSE]
      rownames(mat) <- genes
    }
    .toArrowDF(mat)
  },

  # ---------------------------------------------------------------------------
  # 4. 10x Genomics sparse matrix folder
  # ---------------------------------------------------------------------------
  "10x" = {
    suppressPackageStartupMessages({
      library(Matrix)
      library(Seurat)
    })
    counts <- Read10X(data.dir = inputPath)
    if (normalize10x) {
      obj    <- CreateSeuratObject(counts = counts)
      obj    <- NormalizeData(obj, normalization.method = "LogNormalize",
                              scale.factor = 10000, verbose = FALSE)
      mat    <- GetAssayData(obj, layer = "data")
    } else {
      mat <- counts
    }
    .toArrowDF(mat)
  },

  # ---------------------------------------------------------------------------
  # 5. HDF5 (.h5) — Cell Ranger output
  # ---------------------------------------------------------------------------
  "h5" = {
    suppressPackageStartupMessages({
      library(Seurat)
      library(hdf5r)
    })
    counts <- Read10X_h5(inputPath)
    if (normalizeH5) {
      obj <- CreateSeuratObject(counts = counts)
      obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                           scale.factor = 10000, verbose = FALSE)
      mat <- GetAssayData(obj, layer = "data")
    } else {
      mat <- counts
    }
    .toArrowDF(mat)
  },

  # ---------------------------------------------------------------------------
  # 6. AnnData (.h5ad) — scanpy/Python output
  # ---------------------------------------------------------------------------
  "h5ad" = {
    suppressPackageStartupMessages(library(anndata))
    adata <- read_h5ad(inputPath)
    mat   <- if (is.null(h5adLayer)) {
      t(as.matrix(adata$X))           # X is cells x genes; transpose to genes x cells
    } else {
      t(as.matrix(adata$layers[[h5adLayer]]))
    }
    rownames(mat) <- adata$var_names  # gene names
    colnames(mat) <- adata$obs_names  # cell/sample names
    .toArrowDF(mat)
  },

  stop("Unknown inputFormat: '", inputFormat, "'. ",
       "Choose one of: seurat, sce, tsv, 10x, h5, h5ad")
)

# =============================================================================
# SAVE
# =============================================================================

dir.create(dirname(outputPath), recursive = TRUE, showWarnings = FALSE)
write_feather(df, outputPath)
cat("Saved Arrow file to:", outputPath, "\n")
cat("Dimensions:", nrow(df), "genes x", ncol(df) - 1, "samples/cells\n")
