# Plot weights/confidences
rm(list=ls())
library(ggplot2)

# Custom function to format y-axis labels
custom_labels <- function(x) {
  # Keep 0 as is, scale other values by 1000 and add "K"
  labels <- ifelse(x == 0, "0", paste0(scales::number(x / 1000, accuracy = 0.1), "K"))
  return(labels)
}

histogramConfidencesDir <- function(currNetDirs, breaks) {
    # Ensure breaks length matches currNetDirs length
    if (length(currNetDirs) != length(breaks)) {
        stop("Length of breaks must match length of currNetDirs.")
    }

    for (ix in seq_along(currNetDirs)) {
        currNetDir <- currNetDirs[ix]
        subfolders <- list.dirs(currNetDir, recursive = FALSE)
        target_folders <- subfolders[basename(subfolders) %in% c("TFA", "TFmRNA")]
        
        # Load and plot individually
        for (subfolder in target_folders) {
            # subfolder = target_folders[jx]
            file_path <- file.path(subfolder, "edges_subset.txt")
            if (file.exists(file_path)) {
                data <- read.table(file_path, header = TRUE, sep = "\t")

                # Create individual histogram
                conf <- data.frame(confidence = floor(data$Stability) / max(floor(data$Stability)))
                p <- ggplot((conf), aes(x = confidence)) +
                    geom_histogram(bins = breaks[ix], fill = "blue", color = "black", alpha = 0.9) +
                    labs(title = basename(subfolder), x = "Confidence", y = "Frequency") +
                    theme_bw() +
                    theme(
                        axis.title = element_text(size = 16, color = "black"),
                        axis.text = element_text(size = 14, color = "black"),
                        plot.title = element_text(size = 20,, color = "black", hjust = 0.5) 
                    ) + scale_y_continuous(labels = custom_labels)#+ scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "K"))

                # Save individual plot
                hist_file <- file.path(subfolder, paste0("confidence_distribution_", breaks[ix], ".png"))
                ggsave(hist_file, plot = p, width = 6, height = 5, dpi = 600)
            } else {
                message("File not found: ", file_path)
            }
        }

    }

}


histogramConfidencesData <- function(currNetFiles, breaks) {
  # Ensure breaks length matches currNetFiles length
  if (length(currNetFiles) > length(breaks)) {
    stop("Length of breaks must match or be greater length of currNetFiles.")
  }

  # Loop over the directories in currNetFiles
  for (ix in seq_along(currNetFiles)) {
    currNetFile <- currNetFiles[ix]  # Full path to edges_subset.txt

    # Check if the file exists
    if (file.exists(currNetFile)) {
      data <- read.table(currNetFile, header = TRUE, sep = "\t")

      # Create individual histogram
      conf <- data.frame(confidence = floor(data$Stability) / max(floor(data$Stability)))
      p <- ggplot(conf, aes(x = confidence)) +
        geom_histogram(bins = breaks[ix], fill = "blue", color = "black", alpha = 0.9) +
        labs(title = basename(dirname(currNetFile)), x = "Confidence", y = "Frequency") +
        theme_bw() +
        theme(
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 14, color = "black"),
          plot.title = element_text(size = 20, color = "black", hjust = 0.5)
        ) + scale_y_continuous(labels = custom_labels)  # Custom y-axis labels

      # Save individual plot in the same directory
      hist_file <- file.path(dirname(currNetFile), paste0("confidence_distribution_", breaks[ix], ".png"))
      ggsave(hist_file, plot = p, width = 6, height = 5, dpi = 600)
    } else {
      message("File not found: ", currNetFile)
    }
  }
}

currNetFile <- "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT20/TFmRNA/edges_subset.txt"

conf <- data.frame(confidence = floor(data$Stability))
breaks <- max(conf)
p <- ggplot(conf, aes(x = confidence)) +
  geom_histogram(bins = breaks, fill = "blue", color = "black", alpha = 0.9) +
  labs(title = basename(dirname(currNetFile)), x = "Confidence", y = "Frequency") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 20, color = "black", hjust = 0.5)
  ) 

 hist_file <- file.path(dirname(currNetFile), paste0("testConf_", breaks, ".png"))
 ggsave(hist_file, plot = p, width = 6, height = 5, dpi = 600)

# USAGE


currNetDirs <- c(
    # ---Pseudobulk Inferelator
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63",

    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/Inferelator/ATAC_ChIPprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63",

    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/Inferelator/ATAC_KOprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63",

    # --- single cell inferelator
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p25_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda1p0_220totSS_20tfsPerGene_subsamplePCT10",

    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/SC/lambda0p25_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/SC/lambda1p0_220totSS_20tfsPerGene_subsamplePCT10",

    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/SC/lambda0p25_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/SC/lambda1p0_220totSS_20tfsPerGene_subsamplePCT10",

    # --- single cell downsampled
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27",
    "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13"
)


# Extract subsample fractions
# subsample_fracs <- sapply(currNetDirs, function(x) sub(".*subsampleFrac_([0-9p]+).*", "\\1", x))
# subsample_fracs <- gsub("p", ".", subsample_fracs)  # Convert 'p' to '.' for numeric comparison
# unique_fracs <- unique(subsample_fracs)

breaks = c(rep(150, 9), rep(500, 12))



netFiles =  c(
              "1K" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",
              "10K" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",
              "30K" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13/TFA/edges_subset.txt",
               "77K" = "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10/TFA/edges_subset.txt"
                )
                
histogramConfidencesStacked <- function(currNetFiles, breaks, dirOut, saveName) {

  allData <- data.frame(confidence = numeric(), network = character())

  for (ix in seq_along(currNetFiles)) {
    currNetFile <- currNetFiles[ix] # Full path to edges_subset.txt
    currNetName <- names(currNetFile)
    if (file.exists(currNetFile)) {
      data <- read.table(currNetFile, header = TRUE, sep = "\t")
      # Compute absolute value of signedQuantile and create a temporary data frame.
      conf <- data.frame(confidence = floor(data$Stability) / max(floor(data$Stability)))
      tempDf <- data.frame(confidence = conf,
                          network = currNetName)
      allData <- rbind(allData, tempDf)
    } else {
      message("File not found: ", currNetFile)
    }
  }

  p <- ggplot(allData, aes(x = confidence, color = network)) +
      # geom_freqpoly(bins = breaks[1], size = 1) +
          geom_histogram(bins = breaks[1], alpha = 0.6, position = "identity", color = "black") +
      labs(x = "Confidence", y = "Frequency") +
      theme_bw() +
      theme(
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 14, color = "black")
      )
      histFile <- file.path(dirOut, paste0(saveName, ".pdf"))
      ggsave(histFile, plot = p, width = 6, height = 5, dpi = 600)

}


histogramConfidencesStacked(netFiles, breaks = 300, dirOut = "/data/miraldiNB/Michael/mCD4T_Wayman/Figures", saveName = "confidencesTest")