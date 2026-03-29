# GRN_Tfh10_PR_heatmap_final2.R
# Heatmap of PR results
rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
library(ggplot2)
library(scales)
})
source('/data/miraldiNB/Michael/Scripts/GRN/utils_grn_analysis.R')
source('/data/miraldiNB/wayman/scripts/net_utils.R')
source('/data/miraldiNB/wayman/scripts/clustering_utils.R')

####### INPUTS #######

# Output directory
dir_out <- '/data/miraldiNB/Michael/mCD4T_Wayman/Figures/Tfh10_S4B/2p5PCT'
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# File save name
file_save <- 'Tfh10_S4'

# Gene target file
file_targ <- "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/targRegLists/targetGenes_names.txt"

# TF list
file_tf <- "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/targRegLists/potRegs_names.txt"

# list of networks (in order)
file_grn <- c(

# No Prior : lambdaBias=1 + tfaOpt='_TFmRNA' + ATAC_KOprior
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # NoPrior (TFmRNA + ATAC)

# ATAC Prior
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC)
 "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # TFA 
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv", # +Combine (ATAC)


# ATAC Prior (sc-Inferelator)
#

# ATAC+ChIP Prior
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC + ChIP)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC + ChIP) 
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # TFA 
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC + ChIP)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC + ChIP)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv", # +Combine (ATAC + ChIP)

# ATAC + ChIP Prior (sc- Inferelator)
#

# ATAC+KO Prior
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC + KO)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC + KO) 
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # TFA 
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC + KO)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC + KO)
"/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv", # +Combine (ATAC + KO)

# ATAC + KO Prior (sc- Inferelator)
#

#Cell Oracle
# "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Unprocessed/Combined/Old/max_pCutoff_0.1pct_combined_GRN_3cols.tsv",
"/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Unprocessed/CombinedGRNs/max/max_pCut0.1_meanEdgesPerGene20_combined_GRN.tsv",

# Scenic Plus
"/data/miraldiNB/Michael/mCD4T_Wayman/scenicPlus/eRegulons_direct_filtered.tsv",   # Direct
"/data/miraldiNB/Michael/mCD4T_Wayman/scenicPlus/eRegulons_extended_filtered.tsv", # Extended

#SupirFactor
"/data/miraldiNB/Michael/mCD4T_Wayman/SupirFactor/Hierarchical_20250212_2043/GRN_Hierarchical_long.tsv", # Hierarchical
"/data/miraldiNB/Michael/mCD4T_Wayman/SupirFactor/Shallow_20250330_0212/GRN_Shallow_long.tsv" # Shallow
)

# network names
name_grn <- c(

'No Prior',

'ATAC mRNA +',
'ATAC mRNA ++',
'ATAC TFA',
'ATAC TFA +',
'ATAC TFA ++',
'ATAC Combined +',

# 'sc ATAC mRNA +',
# 'sc ATAC TFA +',
# 'sc ATAC Combined +',

'ChIP mRNA +',
'ChIP mRNA ++',
'ChIP TFA',
'ChIP TFA +',
'ChIP TFA ++',
'ChIP Combined +',

# 'sc ChIP mRNA +',
# 'sc ChIP TFA +',
# 'sc ChIP Combined +',

'KO mRNA +',
'KO mRNA ++',
'KO TFA',
'KO TFA +',
'KO TFA ++',
'KO Combined +',

# 'sc KO mRNA +',
# 'sc KO TFA +',
# 'sc KO Combined +',

'Cell Oracle',
# 'Cell Oracle Mean',

'SCENIC+ Dir',
'SCENIC+ Ext',

'SupirFactor H',
'SupirFactor Sh'

)

# network type
type_grn <- c(

'No Prior',

'ATAC',
'ATAC',
'ATAC',
'ATAC',
'ATAC',
'ATAC',

# 'ATAC',
# 'ATAC',
# 'ATAC',

'ChIP',
'ChIP',
'ChIP',
'ChIP',
'ChIP',
'ChIP',

# 'ChIP',
# 'ChIP',
# 'ChIP',

'KO',
'KO',
'KO',
'KO',
'KO',
'KO',

# 'KO',
# 'KO',
# 'KO',

'Cell Oracle',
# 'Cell Oracle',

'SCENIC+',
'SCENIC+',

'SupirFactor',
'SupirFactor'
)

# Gold standard names
name_gs <- c(
'ChIP',
'KO',
'KC'
)

# List of gold standards
file_gs <- c(
    # '/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/GS/TF_binding/priors/FDR5_Rank50/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv',
    # '/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/GS/RNA/priors/prior_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_sp.tsv',
    # '/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/GS/KC/prior_KC_Thelper_Miraldi2019Th17_Rank100_sp.tsv'
    "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv",
    "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_TF_KO_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_sp.tsv",
    "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_KC_Thelper_Miraldi2019Th17_Rank100_sp.tsv"
)

# gold standard save file name
save_gs <- ''

# params
min_targ_gs <- 20 # gold standard min target cutoff
cutoff_rec <- 0.025 # recall cutoff
range_fc <- c(-2,2) # color scale range for log2fc heatmap
breaks_fc <- c(-2,-1,0,1,2)
heat_h <- 8 # heatmap height
# heat_w <- 15 # heatmap width
heat_w <- 21 # heatmap width
box_w <- 6 # boxplot width

######################

# Output file label
label_pcut <- 100*cutoff_rec
file_save <- if (nzchar(save_gs)) paste0(file_save, '_', save_gs) else file_save


# Load tf list
tf <- readLines(file_tf)

# Calc precision-recall
df_pr <- NULL
for (ix in 1:length(file_grn)){
    
    print(paste0('Process GRN: ',ix,'/',length(file_grn)))
    for (jx in 1:length(file_gs)){
        
        # # Filter gold standard TFs
        gs <- read.delim(file_gs[jx], header=TRUE, sep='\t')
        gs <- subset(gs, TF %in% tf)
        keep_tf <- names(which(table(gs$TF) >= min_targ_gs))
        gs <- subset(gs, TF %in% keep_tf)
        
        # All PR
        curr_tf <- sort(unique(gs$TF))
        pr <- grn_pr(grn=file_grn[ix], gs=gs, gene_tar=file_targ, tf=curr_tf)
        prec <- pr$P; prec <- c(prec, pr$Rand)
        recall <- pr$R; recall <- c(recall, 1)
        
        # Precision at cutoff
        slope_pr <- (prec[which(recall >= cutoff_rec)[1]] - prec[which(recall >= cutoff_rec)[1]-1])/
                        (recall[which(recall >= cutoff_rec)[1]] - recall[which(recall >= cutoff_rec)[1]-1])
        pcut <- prec[which(recall >= cutoff_rec)[1]-1] + 
                        slope_pr*(cutoff_rec-recall[which(recall >= cutoff_rec)[1]-1])
        curr_df <- data.frame(GRN=name_grn[ix], Type=type_grn[ix], GS=name_gs[jx], 
                                TF='All TF', Pcut=pcut, Rand=pr$Rand)
        df_pr <- rbind(df_pr, curr_df)
        
        for (kx in curr_tf){
            
            # pr <- grn_pr(grn=file_grn[ix], gs=file_gs[jx], gene_tar=file_targ, tf=kx)
            pr <- grn_pr(grn=file_grn[ix], gs=gs, gene_tar=file_targ, tf=kx)
            if (!is.null(pr)){
                prec <- pr$P; prec <- c(prec, pr$Rand)
                recall <- pr$R; recall <- c(recall, 1)
                
                # precision at cutoff
                slope_pr <- (prec[which(recall >= cutoff_rec)[1]] - prec[which(recall >= cutoff_rec)[1]-1])/
                                (recall[which(recall >= cutoff_rec)[1]] - recall[which(recall >= cutoff_rec)[1]-1])
                pcut <- prec[which(recall >= cutoff_rec)[1]-1] + 
                                slope_pr*(cutoff_rec-recall[which(recall >= cutoff_rec)[1]-1])
                curr_df <- data.frame(GRN=name_grn[ix], Type=type_grn[ix], GS=name_gs[jx], 
                                        TF=kx, Pcut=pcut, Rand=pr$Rand)
                df_pr <- rbind(df_pr, curr_df)
            }
        }
    }
}
df_pr$PcutLog2FC <- log2(df_pr$Pcut/df_pr$Rand)
df_pr$GS_TF <- paste(df_pr$GS, df_pr$TF, sep='_')

# order gold standards
order_gs <- NULL
for (ix in name_gs){
    order_gs <- c(order_gs, paste(ix,'All TF',sep='_'))
    curr_pr <- subset(df_pr, GS==ix & TF != 'All TF')
    mat <- net_sparse_2_full(curr_pr[,c('GS_TF','GRN','PcutLog2FC')])
    curr_order <- colnames(mat)[cluster_hierarchical(mat)]
    order_gs <- c(order_gs, curr_order)
}
df_pr$GS_TF <- factor(df_pr$GS_TF, levels=order_gs)

# order networks
df_pr$GRN <- factor(df_pr$GRN, levels=rev(name_grn))


# Axis label
label_gs <- order_gs
for(ix in name_gs){
    label_gs <- gsub(paste0(ix,'_'),'',label_gs)
}
face_gs <- rep('plain',length(label_gs))
face_gs[which(label_gs=='All TF')] <- 'bold'

tf_ko <- unique(subset(df_pr, GS=='KO' & TF != 'All TF')$TF)
tf_chip <- unique(subset(df_pr, GS=='ChIP' & TF != 'All TF')$TF)
tf_int <- intersect(tf_ko, tf_chip)
for (ix in tf_int){
    face_gs[which(label_gs==ix)[1]] <- 'italic'
}

# Lines separating gold standards
# vline_sep <- ceiling(cumsum(table(df_pr$GS)[name_gs]/length(name_grn)))+0.5
tf_counts <- tapply(df_pr$TF, df_pr$GS, function(x) length(unique(x)))
vline_sep <- cumsum(tf_counts[name_gs]) + 0.5
vline_sep <- vline_sep[1:(length(vline_sep)-1)]

# lines separating network types
order_type <- rev(unique(type_grn))
hline_sep <- cumsum(table(type_grn)[order_type])+0.5
hline_sep <- hline_sep[1:(length(hline_sep)-1)]

color_pal <-  c(low='dodgerblue3', 'white', high='firebrick3') # c('#2c7bb6','#abd9e9','white','#d7191c','#b10026') 
# Heatmap precision at cutoff
print('Heatmap precision at cutoff')
ggplot(df_pr, aes(x=GS_TF, y=GRN, fill=PcutLog2FC)) + 
            geom_tile() +
            # geom_tile(data = subset(df_pr, !is.na(PcutLog2FC))) +  # Removes NA tiles
            # geom_tile(data = subset(df_pr, is.na(PcutLog2FC)), fill = "lightgrey") +  # Adds NA tiles without borders
            geom_vline(xintercept=vline_sep, lwd=0.8, colour='gray30') +
            geom_hline(yintercept=hline_sep, lwd=0.8, colour='gray30') +
            labs(x=' ', y=' ', fill='Log2 FC') +
            theme(axis.text=element_text(size=12,face='plain')) +
            theme(axis.title=element_text(size=14,face='plain')) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            theme(axis.text.x=element_text(face=face_gs)) +
            # scale_fill_gradient2(low='dodgerblue3', high='firebrick3', limits=range_fc, na.value="lightgrey",oob=squish) + 
            scale_fill_gradientn(
                    # colors=c(low='dodgerblue3', high='firebrick3'),
                    limits=range_fc, 
                    colors=  color_pal, 
                    breaks=breaks_fc, 
                    na.value="lightgrey",oob=squish
                    )+
            scale_x_discrete(labels=label_gs)
file_out <- file.path(dir_out, paste0('heat_',file_save,'_PRcut',label_pcut,'.pdf'))
ggsave(file_out, height=heat_h, width=heat_w)

ggplot(df_pr, aes(x=GS_TF, y=GRN, fill=PcutLog2FC)) + 
            geom_tile() +
            # geom_tile(data = subset(df_pr, !is.na(PcutLog2FC))) +  # Removes NA tiles
            # geom_tile(data = subset(df_pr, is.na(PcutLog2FC)), fill = "lightgrey") +  # Adds NA tiles without borders
            geom_vline(xintercept=vline_sep, lwd=0.8, colour='gray30') +
            geom_hline(yintercept=hline_sep, lwd=0.8, colour='gray30') +
            labs(x=' ', y=' ', fill=' ') +
            theme(axis.text=element_text(size=12,face='plain')) +
            theme(axis.title=element_text(size=14,face='plain')) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            theme(axis.text.x=element_text(face=face_gs)) +
            # scale_fill_gradient2(low='dodgerblue3', high='firebrick3', limits=range_fc, na.value="lightgrey",oob=squish) +
            scale_fill_gradientn(
                    # colors=c(low='dodgerblue3', high='firebrick3'),
                    limits=range_fc, 
                    colors= color_pal, 
                    breaks=breaks_fc, 
                    na.value="lightgrey",oob=squish
                    )+
            scale_x_discrete(labels=label_gs) +
            theme(legend.position='bottom', legend.text=element_text(size=16), legend.key.width=unit(1.5,'cm'))
file_out <- file.path(dir_out, paste0('heat_',file_save,'_PRcut',label_pcut,'_2.pdf'))
ggsave(file_out, height=heat_h, width=heat_w)

# Boxplot log2fc
print('boxplot log2fc')
for (ix in name_gs){
    df_pr_sub <- subset(df_pr, GS==ix & TF != 'All TF')
    ggplot(df_pr_sub, aes(x=GRN, y=PcutLog2FC)) + geom_boxplot() + coord_flip() +
                theme_minimal() +
                labs(x=' ', y='Log2 Fold-Change') +
                theme(axis.text=element_text(size=12,face='plain')) +
                theme(axis.title=element_text(size=14,face='plain')) + 
                theme(axis.line=element_line(size=0.6, colour='grey30', linetype='solid'))
    file_out <- file.path(dir_out, paste0('box_',file_save,'_PRcut',label_pcut,'_',ix,'.pdf'))
    ggsave(file_out, height=heat_h, width=box_w)
}



# set regions where GS is tested on itself to NA
df_pr_grey <- df_pr
df_pr_grey$PcutLog2FC <- ifelse((df_pr$Type == "ChIP" & df_pr$GS != "KO") | 
                            (df_pr$Type == "KO" & df_pr$GS != "ChIP"), 
                            NA, df_pr$PcutLog2FC)

# Heatmap precision at cutoff
print('Heatmap precision at cutoff - greyed out GS tested regions')
ggplot(df_pr_grey, aes(x=GS_TF, y=GRN, fill=PcutLog2FC)) + 
            geom_tile() +
            # geom_tile(data = subset(df_pr, !is.na(PcutLog2FC))) +  # Removes NA tiles
            # geom_tile(data = subset(df_pr, is.na(PcutLog2FC)), fill = "lightgrey") +  # Adds NA tiles without borders
            geom_vline(xintercept=vline_sep, lwd=0.8, colour='gray30') +
            geom_hline(yintercept=hline_sep, lwd=0.8, colour='gray30') +
            labs(x=' ', y=' ', fill='Log2 FC') +
            theme(axis.text=element_text(size=12,face='plain')) +
            theme(axis.title=element_text(size=14,face='plain')) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            theme(axis.text.x=element_text(face=face_gs)) +
            # scale_fill_gradient2(low='dodgerblue3', high='firebrick3', limits=range_fc, na.value="lightgrey",oob=squish) + 
            scale_fill_gradientn(
                    limits=range_fc, 
                    colors=  color_pal, 
                    breaks=breaks_fc, 
                    na.value="lightgrey",oob=squish
                    )+
            scale_x_discrete(labels=label_gs)
file_out <- file.path(dir_out, paste0('heat_',file_save,'_PRcut',label_pcut,'_greyed_out.pdf'))
ggsave(file_out, height=heat_h, width=heat_w)

ggplot(df_pr_grey, aes(x=GS_TF, y=GRN, fill=PcutLog2FC)) + 
            geom_tile() +
            # geom_tile(data = subset(df_pr, !is.na(PcutLog2FC))) +  # Removes NA tiles
            # geom_tile(data = subset(df_pr, is.na(PcutLog2FC)), fill = "lightgrey") +  # Adds NA tiles without borders
            geom_vline(xintercept=vline_sep, lwd=0.8, colour='gray30') +
            geom_hline(yintercept=hline_sep, lwd=0.8, colour='gray30') +
            labs(x=' ', y=' ', fill=' ') +
            theme(axis.text=element_text(size=12,face='plain')) +
            theme(axis.title=element_text(size=14,face='plain')) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            theme(axis.text.x=element_text(face=face_gs)) +
            # scale_fill_gradient2(low='dodgerblue3', high='firebrick3', limits=range_fc, na.value="lightgrey",oob=squish) +
            scale_fill_gradientn(
                    # colors=c(low='dodgerblue3', high='firebrick3'),
                    limits=range_fc, 
                    colors=  color_pal, 
                    breaks=breaks_fc, 
                    na.value="lightgrey",oob=squish
                    )+
            scale_x_discrete(labels=label_gs) +
            theme(legend.position='bottom', legend.text=element_text(size=16), legend.key.width=unit(1.5,'cm'))
file_out <- file.path(dir_out, paste0('heat_',file_save,'_PRcut',label_pcut,'_greyed_out_2.pdf'))
ggsave(file_out, height=heat_h, width=heat_w)

# Boxplot log2fc
print('boxplot log2fc')
for (ix in name_gs){
    df_pr_sub <- subset(df_pr_grey, GS==ix & TF != 'All TF')
    ggplot(df_pr_sub, aes(x=GRN, y=PcutLog2FC)) + geom_boxplot() + coord_flip() +
                theme_minimal() +
                labs(x=' ', y='Log2 Fold-Change') +
                theme(axis.text=element_text(size=12,face='plain')) +
                theme(axis.title=element_text(size=14,face='plain')) + 
                theme(axis.line=element_line(size=0.6, colour='grey30', linetype='solid'))
    file_out <- file.path(dir_out, paste0('box_',file_save,'_PRcut',label_pcut,'_',ix,'_greyed_out.pdf'))
    ggsave(file_out, height=heat_h, width=box_w)
}                 
print('DONE')
