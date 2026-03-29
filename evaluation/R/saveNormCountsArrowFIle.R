rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(arrow)  
})

# make sure you load or have the latest version of curl installed.

# make sure you load or have the latest version of curl installed.

rds_path <- "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/annotation_scrna_final/obj_Tfh10_RNA_annotated.rds"
obj <- readRDS(rds_path)

# Process File for scGRN using Inferelator.
obj <- NormalizeData(obj)
norm_counts <- as.matrix(GetAssayData(obj, layer='data'))
# write_parquet(as.data.frame(norm_counts), "/data/miraldiNB/Michael/GRN_Benchmark/Data/Tfh10_scRNA_logNorm_Counts.parquet")

norm_counts <- as.data.frame(norm_counts)
norm_counts$Genes <- rownames(norm_counts)
norm_counts <- norm_counts[, c("Genes", setdiff(colnames(norm_counts), "Genes"))] # Reorder the columns to make 'Genes' the first column

# Save as an arrow file
write_feather(norm_counts, "/data/miraldiNB/Michael/GRN_Benchmark/Data/Tfh10_scRNA_logNorm_Counts.arrow")

# write_feather(norm_counts, "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/scNormCounts.arrow")

# feather_file <- read_feather("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Pseudobulk/RNA2/SC_counts.feather")
# write_ipc_file(feather_file, "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Pseudobulk/RNA2/SC_counts.arrow")


rds_path <- "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/ObjFiltered.rds"