# =============================================================================
#  validateBEBICLambdas.R — Compare per-gene EBIC-optimal lambda distributions
#                           across bEBIC runs (auto vs 40-step vs 100-step grid)
#
#  Purpose:
#    Confirms that user-supplied lambda grids (40 or 100 steps) resolve the same
#    lambda region as GLMNet's automatic path. If the per-gene median lambda ECDFs
#    overlap, the grids are equivalent and network differences are not an artefact
#    of grid resolution.
#
#  Input:  bebic_lambda_summary.tsv (gene | median_lambda | std_lambda)
#          Written by bebicFinalFit!(... outputDir=...) — re-run with outputDir
#          set for any run that is missing this file (see NOTE below).
#
#  Outputs (saved to dirOut):
#    lambdaECDF_medianLambda_<subfolder>.pdf  — ECDF of per-gene median lambda
#    lambdaDensity_stdLambda_<subfolder>.pdf  — density of per-gene lambda std
#    lambdaSummary_stats_<subfolder>.tsv      — numeric summary table
#
#  NOTE: To generate bebic_lambda_summary.tsv for a run that is missing it,
#        re-run bebicFinalFit! with outputDir set:
#          bebicFinalFit!(grnData, buildGrn; outputDir = "/path/to/run/TFA")
#        This is cheap — it does NOT re-run the subsampling stage.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})


# =============================================================================
# USER INPUTS
# =============================================================================

# Named list: label => path to bebic_lambda_summary.tsv
# Add or remove runs freely — missing files are skipped with a warning.
lambdaFiles <- list(
  "bEBIC_autoLambda"    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max/TFA/bebic_lambda_summary.tsv",
  "bEBIC_40lambda"      = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_lambdaProvided/TFA/bebic_lambda_summary.tsv",
  "bEBIC_100lambda"     = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_100lambdaProvided/TFA/bebic_lambda_summary.tsv",
  "bEBIC_noZ_40lambda"  = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_noZscoreTFA_noZscoreLASSO_lambdaProvided/TFA/bebic_lambda_summary.tsv",
  "bEBIC_noZ_100lambda" = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_noZscoreTFA_noZscoreLASSO_100lambdaProvided/TFA/bebic_lambda_summary.tsv"
)

# Subfolder label used in output filenames (e.g. "TFA" or "TFmRNA")
subfolderLabel <- "TFA"

# Output directory
dirOut <- "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/netEval/lambdaValidation"
dir.create(dirOut, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# LOAD — skip missing files with an informative message
# =============================================================================

lambdaDfList <- lapply(names(lambdaFiles), function(nm) {
  path <- lambdaFiles[[nm]]
  if (!file.exists(path)) {
    message(sprintf(
      "MISSING: %s\n  -> %s\n  -> Re-run bebicFinalFit!(grnData, buildGrn; outputDir=\"%s\") to generate it.",
      nm, path, dirname(path)
    ))
    return(NULL)
  }
  df      <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$run  <- nm
  df
})

lambdaDf <- do.call(rbind, Filter(Negate(is.null), lambdaDfList))

if (is.null(lambdaDf) || nrow(lambdaDf) == 0) {
  stop("No lambda summary files found. See messages above for how to generate them.")
}

runsLoaded  <- unique(lambdaDf$run)
runsMissing <- setdiff(names(lambdaFiles), runsLoaded)

message(sprintf("\nLoaded %d run(s): %s", length(runsLoaded), paste(runsLoaded, collapse = ", ")))
if (length(runsMissing) > 0) {
  message(sprintf("Skipped %d run(s) (missing file): %s\n",
                  length(runsMissing), paste(runsMissing, collapse = ", ")))
}

# Fix run order for legend
lambdaDf$run <- factor(lambdaDf$run, levels = names(lambdaFiles)[names(lambdaFiles) %in% runsLoaded])

# Colour palette — one colour per run
runColors <- setNames(
  scales::hue_pal()(length(runsLoaded)),
  runsLoaded
)


# =============================================================================
# SHARED THEME
# =============================================================================

gridTheme <- theme(
  panel.grid.major = element_line(color = "grey75",  linewidth = 0.3),
  panel.grid.minor = element_line(color = "#e8e8e8", linewidth = 0.1),
  panel.border     = element_rect(color = "grey50",  linewidth = 0.3, fill = NA),
  axis.ticks       = element_line(linewidth = 0.15, color = "grey60"),
  axis.ticks.length = unit(2, "pt")
)


# =============================================================================
# PLOT 1 — ECDF of per-gene median lambda (log scale)
#
# Story: if the ECDFs from auto / 40-step / 100-step overlap, the user-supplied
# grids are resolving the same lambda region as GLMNet's automatic path.
# A shifted curve means the grid boundaries were too narrow (confirm with
# the boundary hit warnings from bebicEstimateLambdas!).
# =============================================================================

p_ecdf <- ggplot(lambdaDf, aes(x = median_lambda, color = run)) +
  stat_ecdf(linewidth = 0.7) +
  scale_x_log10(
    labels       = scales::label_number(accuracy = 0.001),
    breaks       = 10^seq(-4, 1),
    minor_breaks = as.vector(outer(1:9, 10^seq(-4, 1)))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent
  ) +
  scale_color_manual(values = runColors) +
  annotation_logticks(sides     = "b",
                      linewidth = 0.2,
                      color     = "grey60",
                      short     = unit(1.5, "pt"),
                      mid       = unit(2.5, "pt"),
                      long      = unit(3.5, "pt")) +
  labs(x     = "Per-gene median lambda (log scale)",
       y     = "Cumulative proportion of genes",
       title = paste0("bEBIC lambda ECDF — ", subfolderLabel)) +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text       = element_text(size = 7),
        axis.title      = element_text(size = 9),
        legend.position = "bottom",
        legend.text     = element_text(size = 7),
        legend.title    = element_blank()) +
  gridTheme

ggsave(file.path(dirOut, paste0("lambdaECDF_medianLambda_", subfolderLabel, ".pdf")),
       plot = p_ecdf, width = 5, height = 4, dpi = 600)
message("Saved: lambdaECDF_medianLambda_", subfolderLabel, ".pdf")


# =============================================================================
# PLOT 2 — Density of per-gene lambda std (log scale)
#
# Story: high std means the EBIC-optimal lambda varies a lot across subsamples
# for that gene — unstable selection. Comparing std distributions tells you
# whether one grid setting leads to more variable lambda selection than another.
# Expected: z-scored runs have lower std (lambda scale is normalised).
# =============================================================================

# Drop genes where std == 0 (all subsamples chose the same lambda — can happen
# at the edge of the grid; also causes -Inf on log scale)
lambdaDf_std <- lambdaDf[lambdaDf$std_lambda > 0, ]

p_std <- ggplot(lambdaDf_std, aes(x = std_lambda, color = run, fill = run)) +
  geom_density(alpha = 0.15, linewidth = 0.6) +
  scale_x_log10(
    labels       = scales::label_number(accuracy = 0.001),
    breaks       = 10^seq(-4, 1),
    minor_breaks = as.vector(outer(1:9, 10^seq(-4, 1)))
  ) +
  scale_color_manual(values = runColors) +
  scale_fill_manual(values  = runColors) +
  annotation_logticks(sides     = "b",
                      linewidth = 0.2,
                      color     = "grey60",
                      short     = unit(1.5, "pt"),
                      mid       = unit(2.5, "pt"),
                      long      = unit(3.5, "pt")) +
  labs(x     = "Per-gene lambda std across subsamples (log scale)",
       y     = "Density",
       title = paste0("bEBIC lambda variability — ", subfolderLabel)) +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text       = element_text(size = 7),
        axis.title      = element_text(size = 9),
        legend.position = "bottom",
        legend.text     = element_text(size = 7),
        legend.title    = element_blank()) +
  gridTheme

ggsave(file.path(dirOut, paste0("lambdaDensity_stdLambda_", subfolderLabel, ".pdf")),
       plot = p_std, width = 5, height = 4, dpi = 600)
message("Saved: lambdaDensity_stdLambda_", subfolderLabel, ".pdf")


# =============================================================================
# NUMERIC SUMMARY TABLE
#
# Median, IQR, and min/max of per-gene median_lambda per run.
# Quick sanity check that the lambda scales are in the same ballpark.
# =============================================================================

summaryStats <- lambdaDf %>%
  group_by(run) %>%
  summarise(
    n_genes       = n(),
    lambda_min    = min(median_lambda),
    lambda_q25    = quantile(median_lambda, 0.25),
    lambda_median = median(median_lambda),
    lambda_q75    = quantile(median_lambda, 0.75),
    lambda_max    = max(median_lambda),
    std_median    = median(std_lambda),
    .groups       = "drop"
  )

print(summaryStats, n = Inf)
write.table(summaryStats,
            file.path(dirOut, paste0("lambdaSummary_stats_", subfolderLabel, ".tsv")),
            quote = FALSE, row.names = FALSE, sep = "\t")
message("Saved: lambdaSummary_stats_", subfolderLabel, ".tsv")
