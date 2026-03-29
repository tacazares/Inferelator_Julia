rm(list = ls())
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(dplyr)
    library(tidyverse)
    library(circlize)
    library(RColorBrewer) 
    library(gridtext)
    library(readxl)
    library(reshape2)
    library(Seurat)
    library(ggplot2)
})

source("/data/miraldiNB/Michael/Scripts/standardize_normalize_pseudobulk.R")
source("/data/miraldiNB/Michael/Scripts/GSEA/GSEA_utils.R")

dirOut <- "/data/miraldiNB/Michael/mCD4T_Wayman/Figures/Fig4/Test"
dir.create(dirOut, showWarnings = F, recursive = T)


# network <- "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.tsv"
tfa_file <- "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/TFA.txt"
file_k_clust_across <- NULL #file.path(dir_out, 'k_clust_across.rds')
file_k_clust_within <- NULL #file.path(dir_out, 'k_clust_within.rds') 

order_celltype <- c('Tfh10','Tfh_Int','Tfh',
                        'Tfr','cTreg','eTreg','rTreg','Treg_Rorc',
                        'Th17','Th1','CTL_Prdm1','CTL_Bcl6',
                        #'TM_Act','TM_ISG',
                        'TEM','TCM')
order_age <- c('Young','Old')
order_rep <- c('R1','R2','R3','R4')

tfa <- read.table(tfa_file)
meta_data <- NULL
for (ix in order_celltype){
    for (jx in order_age){
        curr_sample <- paste(order_rep, jx, ix, sep='_')
        curr_meta <- data.frame(row.names=curr_sample, CellType=ix, Age=jx, Rep=order_rep)
        meta_data <- rbind(meta_data, curr_meta)
    }
}

# diff <- setdiff(rownames(meta_data), colnames(tfa))
meta_data <- meta_data[intersect(colnames(tfa), rownames(meta_data)), ]

# rearrange tfa matrix
tfa <- tfa[ ,rownames(meta_data)]
# Z-scoring or Standardization
z_tfa <- standardizeAndNormalizeCounts(tfa, meta_data = meta_data, celltype_var = "CellType", epsilon = NULL)
tfa_z_across = z_tfa$z_across
tfa_z_within = z_tfa$z_within

# Create Annotation Color
getPalette = colorRampPalette(brewer.pal(13, "Set1"))
celltype_colors = getPalette(length(unique(meta_data[,1])))
celltype_colors <- setNames(celltype_colors, c(unique(meta_data[,1])))
treatment_colors <- c('Young'='grey66', 'Old'='grey33')

meta_data[,2] <- factor()
#Create annotation column
colAnn_top <- HeatmapAnnotation(
                    `Cell Type` = meta_data[,1],
                    `Age` = meta_data[,2],
                    col = list('Cell Type' = celltype_colors,
                              'Age' = treatment_colors),
                              show_annotation_name = FALSE
                        )   

across_center <- 6
# within_center <- 4
gc()
if (is.null(file_k_clust_across)){
  k_clust_across <- kmeans(tfa_z_across, centers=across_center,nstart=20, iter.max = 50)
}else{
  k_clust_across <- readRDS(file_k_clust_across)
}

# gc()
# if (is.null(file_k_clust_within)){
#   k_clust_within <- kmeans(tfa_z_within, centers=within_center, nstart=20, iter.max = 50)
# }else{
#   k_clust_within <- readRDS(file_k_clust_within)
# }

# # Cluster
k_clust_across <- kmeans(tfa_z_across, centers=6,nstart=20, iter.max = 50)
split_by <- factor(k_clust_across$cluster, levels = c(1:6))

heat_col <- colorRamp2(c(-2, 0, 2), c('dodgerblue3','white','red')) ## Colors for heatmap
pdf(file.path(dirOut,  "htmap_zscore_across.pdf"), width = 4, height = 6, compress = T)
ht1 <- Heatmap(tfa_z_across,
            name='Z-score',
            show_row_names = F,
            row_split=split_by,
            row_gap = unit(0, "mm"),
            column_gap = unit(0, "mm"), 
            border = TRUE,
            row_title = NULL,
            row_dend_reorder = F,
            use_raster = F,
            col = heat_col,
            # clustering settings
            cluster_rows = TRUE,                    # allow hierarchical clustering
            cluster_row_slices = TRUE,              # reorder within each k-means cluster
            cluster_columns = FALSE,
            cluster_column_slices = F,
            # column_order = colnames(tfa_z_across),

            # row/column dendrograms
            show_row_dend = FALSE,                   # show dendrogram within slices
            show_column_dend = FALSE,

            show_column_names = FALSE,
            column_names_side = 'top',
            # column_names_rot = 45,
            column_title = NULL,
            column_split = meta_data[,1],
            top_annotation = colAnn_top
              #bottom_annotation = colAnn_bottom
              )
draw(ht1)
dev.off()


# cluster_annot <- data.frame(across=k_clust_across$cluster,within=k_clust_within$cluster)
# cluster_annot$cluster <- paste0(cluster_annot$across,cluster_annot$within)
# cluster_annot <- cluster_annot[order(cluster_annot$cluster), ]
# within_clust_order <- c(1,2,3,4)
# across_clust_order <- c(2,3,5,4,6,1)
# cluster_annot$within <- factor(cluster_annot$within, levels=within_clust_order, ordered=T)
# cluster_annot$across <- factor(cluster_annot$across, levels=across_clust_order, ordered=T) 


# annot_cols_across <- c('1'='#cc0001','2'='#fb940b','3'='#ffff01','4'='#01cc00','5'='#2085ec','6'='#fe98bf','7'='#762ca7','8'='#ad7a5b', 
#                     '9'='grey50', '10'='turquoise4')

# annot_cols_within <- c('1'='#badf55','2'='#35b1c9','3'='#b06dad','4'="#14A76C", '5' = 'grey90', '6' = '#e96060')

# #reorder clusters
# split <- factor(cluster_annot$cluster, levels=c(
#                                                 21,22,23,24,
#                                                 31,32,33,34,
#                                                 51,52,53,54,
#                                                 41,42,43,44,
#                                                 61,62,63,64,
#                                                 11,12,13,14
#                                                  ))

# cluster_annot$cluster <- factor(cluster_annot$cluster, levels=levels(split))
# cluster_annot <- cluster_annot[order(cluster_annot$cluster),]
# names(split) <- rownames(cluster_annot)

# rowAnn_across_left1 <- HeatmapAnnotation(`Across` = cluster_annot[,'across'], #df = cluster_annot[,'across'],
#                                     col = list('Across'= annot_cols_across),
#                                     which = 'row',
#                                     simple_anno_size = unit(3, 'mm'),
#                                     #annotation_width = unit(c(1, 4), 'cm'),
#                                     #gap = unit(0, 'mm'),
#                                     show_annotation_name = F)


# pdf(file.path(dirOut,  "htmap_zscore_across1.pdf"), width = 4, height = 6, compress = T)
# ht1 <- Heatmap(tfa_z_across[rownames(cluster_annot),],
#               name='Z-score',
#               show_row_names = F,
#               row_split=cluster_annot$across,
#               show_column_dend =F,
#               row_gap = unit(0, "mm"),
#               column_gap = unit(0, "mm"), 
#               border = TRUE,
#               row_title = NULL,
#               row_dend_reorder = F,
#               use_raster = F,
#               col = heat_col,
#               column_names_gp = gpar(fontsize = 8, fontface = 2),
#               show_row_dend = FALSE,
#               cluster_rows = TRUE,
#               cluster_columns = F,
#               cluster_row_slices = F,
#               cluster_column_slices = F,
#               column_order = colnames(tfa_z_across),
#               show_column_names = FALSE,
#               column_names_side = 'top',
#             #   column_labels = gt_render(column_label),
#               column_names_rot = 45,
#               column_title = NULL,
#               column_split = meta_data[,1],
#               left_annotation = rowAnn_across_left1,
#             #   right_annotation = rowAnn_across_right1,
#               row_order = rownames(cluster_annot),
#               top_annotation = colAnn_top
#               #bottom_annotation = colAnn_bottom
#               )



# PART C: TFA Visualization between data representation
tfaFiles <- list(
            "PB" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/Bulk/ATAC/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/TFA.txt",
            "SC" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/SC/ATAC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5/Combined/TFA.txt",
            "MC2" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/metaCells/MC2/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63_logNorm/Combined/TFA.txt",
            "SEA" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/metaCells/SEACells/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/Combined/TFA.txt"
)

tfaMatList <- list()
for (typeName in names(tfaFiles)) {
  filePath <- tfaFiles[[typeName]]
  tfaMatList[[typeName]] <- read.table(filePath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

numTFA <- length(tfaMatList)
tfaNames <- names(tfaMatList)



# -----------------------------
# Feature Plot
# -----------------------------
objPath <- "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/annotation_scrna_final/obj_Tfh10_RNA_annotated.rds"
tfa_file1 <- "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/noMergedTF/SC/ATAC/geneLambda0p5_220totSS_20tfsPerGene_subsamplePCT5/Combined/TFA.txt"
obj <- readRDS(objPath)
tfa_mat1 <- read.table(tfa_file1, header = T)
lineage_TFs <- c("Bcl6","Maf","Batf","Tox","Ascl2",
                 "Foxp3","Ikzf2","Rorc",
                 "Rorc","Stat3",
                 "Tbx21","Stat4","Eomes","Prdm1",
                 "Tcf7","Klf2","Lef1")

# Filter lineage TFs present in your TFA matrix
lineage_TFs_present <- intersect(lineage_TFs, rownames(tfa_mat))

# Add each TFA as metadata
for(tf in lineage_TFs_present){
  obj[[tf]] <- tfa_mat[tf, colnames(obj)]
}

# Example: FeaturePlot for all lineage TFs
FeaturePlot(seurat_obj,
            features = lineage_TFs_present,
            cols = c("lightgrey","red"),
            reduction = "umap")  # or "tsne" if you use tSNE

