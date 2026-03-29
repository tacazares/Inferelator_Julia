library(biomaRt)
library(stringr)
library(DESeq2)


counts <- read.table("GSE271788_dedup_counts.txt", header = T)
rownames(counts) <- counts[,1]
counts <- counts[,7:length(counts)]

new_names <- str_extract(colnames(counts), "Donor_\\d+_[A-Za-z0-9]+")
new_names <- sapply(strsplit(new_names, "_"), function(x) paste(x[3], x[1], x[2], sep = "_"))
colnames(counts) <- new_names

ensembl_ids_clean <- sub("\\..*", "", rownames(counts))
gtf <- import("gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz")
# Keep only gene entries
genes <- gtf[gtf$type == "gene"]

# Build mapping table
mapping <- data.frame(
  ensembl_gene_id = gsub("\\..*", "", genes$gene_id),  # remove version
  gene_name = genes$gene_name
) %>%
  distinct()
symbol_map <- mapping$gene_name[match(gene_ids, mapping$ensembl_gene_id)]
symbol_map[is.na(symbol_map)] <- gene_ids[is.na(symbol_map)]


gene_ids <- gsub("\\..*", "", rownames(counts))  # remove versions
symbol_map <- mapping$gene_name[match(gene_ids, mapping$ensembl_gene_id)]
rownames(counts) <- ifelse(is.na(symbol_map), gene_ids, symbol_map)
counts <- counts[!grepl("^ENSG", rownames(counts)), ]

# Build metadata dataframe
colnames(counts) <- make.unique(colnames(counts))

clean_names <- sub("\\.\\d+$", "", colnames(counts))

# Split by underscore
split_info <- do.call(rbind, strsplit(clean_names, "_"))

# Build metadata dataframe
metadata <- data.frame(
  TF = split_info[, 1],
  Donor = paste(split_info[, 2], split_info[, 3], sep = "_"),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Donor + TF
)
dds <- DESeq(dds)

#cor_matrix <- dcast(network, Gene ~ TF, value.var = "signedQuantile", fill = 0)


### Plot TFA changes
library(tidyverse)
bulk_cor <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/FullPseudobulk/FullPseudobulk/lambda0p5_200totSS_20tfsPerGene_subsamplePCT63/Combined/TFA_cor.txt")
bulk_quant <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/FullPseudobulk/FullPseudobulk/lambda0p5_200totSS_20tfsPerGene_subsamplePCT63/Combined/TFA_quant.txt")
bulk_sign <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/FullPseudobulk/FullPseudobulk/lambda0p5_200totSS_20tfsPerGene_subsamplePCT63/Combined/TFA_sign.txt")
sc_cor <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/scSubsampleFraction/lambda0p5_200totSS_20tfsPerGene_subsamplePCT10/Combined/TFA_cor.txt")
sc_quant <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/scSubsampleFraction/lambda0p5_200totSS_20tfsPerGene_subsamplePCT10/Combined/TFA_quant.txt")
sc_sign <- read.table("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/scSubsampleFraction/lambda0p5_200totSS_20tfsPerGene_subsamplePCT10/Combined/TFA_binary.txt")
tfa_matrix <- list(bulk_cor, bulk_quant, bulk_sign, sc_cor, sc_quant, sc_sign)
names(tfa_matrix) <- c("bulk_cor","bulk_quant","bulk_sign","sc_cor","sc_quant","sc_sign")

tf_list <- c(
  "AIRE", "BACH2", "BCL11B", "BPTF", "CLOCK", "EPAS1", "ETS1", "FOXK1", "FOXP1", "FOXP3", 
  "GATA3", "GFI1", "HIVEP2", "IKZF1", "IRF1", "IRF2", "IRF4", "IRF7", "IRF9", "KLF2", 
  "KMT2A", "MBD2", "MYB", "NFE2L2", "NFAT5", "NFKB1", "NFKB2", "RELA", 
  "RELB", "REL", "RFX5", "RORC", "SETDB1", "SREBF1", "STAT1", "STAT2", "STAT3", "STAT5A", 
  "STAT5B", "TBX21", "TCF3", "TP53", "YBX1", "YBX3", "YY1", "ZBTB14", "ZFP3", "ZKSCAN1", 
  "ZNF329", "ZNF341", "ZNF791"
)
tf_list <- tf_list[which(tf_list %in% rownames(tfa_matrix[[1]]))]
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


analyse_one_matrix <- function(tfa_matrix,
                               tf_list,
                               pdf_name,
                               sig_cutoffs = c(`***` = 0.001,
                                               `**`  = 0.01,
                                               `*`   = 0.05)) {

  p_to_symbol <- function(p) {
    stars <- names(sig_cutoffs)[which(p < sig_cutoffs)]
    if (length(stars) == 0) "ns" else stars[1]
  }

  res <- purrr::map_dfr(tf_list, function(tf) {

    #--- pull the row safely ---------------------------------------------------
    if (!tf %in% rownames(tfa_matrix)) {
      warning("TF ", tf, " not present in matrix; skipping.")
      return(tibble::tibble(TF = tf, p.value = NA_real_,
                            signif = "NA", direction = "NA", plot = list(NULL)))
    }

    tf_activity <- tfa_matrix[tf, ]           # numeric vector
    tf_df <- tibble::tibble(Sample   = names(tf_activity),
                            Activity = as.numeric(tf_activity)) %>%
      dplyr::mutate(
        Knockout_TF = sub("_Donor_.*", "", Sample),
        Donor       = sub(".*_Donor_", "", Sample) %>% sub("\\..*", "", .)
      )

    #--- average per donor -----------------------------------------------------
    controls  <- tf_df %>% dplyr::filter(Knockout_TF == "AAVS1") %>%
      dplyr::group_by(Donor) %>%
      dplyr::summarise(Control = mean(Activity), .groups = "drop")

    knockouts <- tf_df %>% dplyr::filter(Knockout_TF == tf) %>%
      dplyr::group_by(Donor) %>%
      dplyr::summarise(Knockout = mean(Activity), .groups = "drop")

    paired <- dplyr::inner_join(controls, knockouts, by = "Donor")

    if (nrow(paired) < 2) {
      warning("TF ", tf, ": fewer than 2 paired donors – skipped")
      return(tibble::tibble(TF        = tf, p.value = NA_real_,
                            signif    = "NA", direction = "NA", plot = list(NULL)))
    }

    #--- stats -----------------------------------------------------------------
    t_res      <- stats::t.test(paired$Knockout, paired$Control, paired = TRUE)
    p_val      <- t_res$p.value
    direction  <- ifelse(mean(paired$Knockout - paired$Control) > 0, "Up", "Down")
    signif_sym <- p_to_symbol(p_val)

    #--- plot ------------------------------------------------------------------
    long <- paired %>%
      tidyr::pivot_longer(Control:Knockout,
                          names_to = "Condition",
                          values_to = "Activity")

    g <- ggplot2::ggplot(long, ggplot2::aes(Condition, Activity, group = Donor)) +
      ggplot2::geom_line(ggplot2::aes(color = Donor), linewidth = 1.1, alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Donor), size = 3) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        title     = paste(tf, "Activity Before and After Knockout"),
        subtitle  = paste0("p = ", signif(p_val, 3),
                           " (", signif_sym, ", ", direction, ")"),
        x = "Condition", y = "Estimated TF Activity", color = "Donor") +
      ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))

    tibble::tibble(TF = tf,
                   p.value = p_val,
                   signif  = signif_sym,
                   direction = direction,
                   plot = list(g))
  })

  #--- write one PDF with all TFs ---------------------------------------------
  grDevices::pdf(pdf_name, width = 8, height = 10)
  purrr::walk(res$plot[!purrr::map_lgl(res$plot, is.null)], print)
  grDevices::dev.off()

  res            # keep the plot column this time
}

################################################################
## 2.  MANY-MATRIX WRAPPER (robust, no imap_dfr required)     ##
################################################################
analyse_many_matrices <- function(tfa_list,   # named list or named file vector
                                  tf_list,
                                  out_dir = ".") {

  # read .rds files if character vector of paths was supplied
  if (is.character(tfa_list) && !is.matrix(tfa_list[[1]])) {
    tfa_list <- lapply(tfa_list, readRDS)
    names(tfa_list) <- basename(names(tfa_list))    # keep vector names
  }

  if (is.null(names(tfa_list)) || any(names(tfa_list) == ""))
    stop("tfa_list must be a *named* list or vector so datasets can be labelled.")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  results <- vector("list", length(tfa_list))
  out_i   <- 1L

  for (mat_name in names(tfa_list)) {

    mat <- tfa_list[[mat_name]]
    if (!(is.matrix(mat) || is.data.frame(mat))) {
      warning("Skipping '", mat_name, "': not a matrix/data.frame")
      next
    }

    pdf_file <- file.path(out_dir,
                          paste0("TF_activity_changes_", mat_name, ".pdf"))

    message("Processing ", mat_name, " …")
    stats_one <- tryCatch(
      analyse_one_matrix(mat, tf_list, pdf_name = pdf_file),
      error = function(e) {
        warning("Failed on ", mat_name, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(stats_one)) next         # skip failures

    stats_one <- tibble::as_tibble(stats_one)        # force tibble
    stats_one <- dplyr::mutate(stats_one, dataset = mat_name, .before = 1)

    results[[out_i]] <- stats_one
    out_i <- out_i + 1L
  }

  dplyr::bind_rows(results[seq_len(out_i - 1L)])
}


all_stats <- analyse_many_matrices(tfa_matrix,
                                   tf_list,
                                   out_dir = "plots")   # PDFs in ./plots

##  View or export the combined statistics:
print(all_stats)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(tidyr)

pct_expressing <- function(obj,
                           genes,
                           assay = "RNA",
                           slot  = "data") {

  expr <- Seurat::GetAssayData(obj[[assay]], slot = slot)   # genes × cells
  vec  <- vapply(genes, function(g) {
             if (!g %in% rownames(expr)) return(NA_real_)
             mean(expr[g, ] > 0) * 100                      # % of cells
           }, numeric(1))
  names(vec) <- genes
  vec
}

DefaultAssay(obj) <- "RNA"               # if not already set
pct_vec <- pct_expressing(obj, tf_list) 

score_mat <- all_stats %>%
  mutate(score = ifelse(
    is.na(p.value),
    NA_real_,
    -log10(p.value) * ifelse(direction == "Up",  1, -1))
  ) %>%
  select(TF, dataset, score) %>%
  pivot_wider(names_from = dataset, values_from = score) %>%
  as.data.frame()

rownames(score_mat) <- score_mat$TF
score_mat <- as.matrix(score_mat[, -1, drop = FALSE])

star_mat <- all_stats %>%
  mutate(stars = case_when(
    is.na(p.value)      ~ "",
    p.value < 0.001     ~ "***",
    p.value < 0.01      ~ "**",
    p.value < 0.10      ~ "*",
    TRUE                ~ ""
  )) %>%
  select(TF, dataset, stars) %>%
  pivot_wider(names_from = dataset, values_from = stars) %>%
  as.data.frame()

rownames(star_mat) <- star_mat$TF
star_mat <- as.matrix(star_mat[, -1, drop = FALSE])

deg_counts <- setNames(rep(0, nrow(score_mat)), rownames(score_mat))
deg_counts[names(target_number)] <- target_number      # overwrite where present

# ---- LEFT annotation: % of cells expressing ----------------------------------
row_anno_left <- rowAnnotation(
  `% Cells\nExpressing` = anno_barplot(
    pct_vec[rownames(score_mat)],                 # ensure correct order
    gp          = gpar(fill = "steelblue", col = NA),
    bar_width   = 0.85,
    border      = FALSE,
    axis_param  = list(at = c(0, 50, 100), gp = gpar(fontsize = 7))
  ),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  width = unit(2.4, "cm")
)

# ---- RIGHT annotation: # of DE genes -----------------------------------------
row_anno_right <- rowAnnotation(
  `# DE genes` = anno_barplot(
    deg_counts,
    gp          = gpar(fill = "grey40", col = NA),
    bar_width   = 0.85,
    border      = FALSE,
    axis_param  = list(at = c(0, max(deg_counts)), gp = gpar(fontsize = 7))
  ),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  width = unit(2.3, "cm")
)

max_abs  <- max(abs(score_mat), na.rm = TRUE)
col_fun  <- circlize::colorRamp2(c(-max_abs, 0, max_abs),
                                 c("navy", "white", "firebrick"))
text_col <- function(fill) {
  rgb <- col2rgb(fill) / 255
  ifelse(0.299*rgb[1] + 0.587*rgb[2] + 0.114*rgb[3] < 0.5, "white", "black")
}

row_labels <- rowAnnotation(
  TF = anno_text(
    rownames(score_mat),
    gp        = gpar(fontsize = 9),
    just      = "left",
    location  = 0.5        # centred vertically in each cell
  ),
  width = unit(2.2, "cm"), # reserve space for the longest name
  show_annotation_name = FALSE
)

###############################################################################
## 2.  Core heat-map (row names turned OFF)  ##################################
###############################################################################
ht_core <- Heatmap(
  score_mat,
  name               = "-log10(p)",
  col                = col_fun,
  na_col             = "grey90",
  border             = TRUE,
  row_km             = 3,

  show_row_names     = FALSE,      # <── row names handled by row_labels
  column_names_gp    = gpar(fontsize = 9),

  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- star_mat[i, j]
    if (nzchar(lab)) {
      grid.text(
        lab, x, y,
        gp = gpar(fontsize = 8,
                  fontface = "bold",
                  col = text_col(fill))
      )
    }
  }
)

###############################################################################
## 3.  Assemble  LEFT  +  CORE  +  LABELS  +  RIGHT  ##########################
###############################################################################
ht_full <- row_anno_left + ht_core + row_labels + row_anno_right
#            %-expressing      matrix     TF names      #DE genes

###############################################################################
## 4.  Draw / save  ###########################################################
###############################################################################
pdf("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/subNetworks/scSubsampleFraction/lambda0p5_200totSS_20tfsPerGene_subsamplePCT10/Combined/plots//TF_heatmap_with_two_bars_and_labels.pdf", width = 7.5, height = 9)
draw(ht_full,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")
dev.off()