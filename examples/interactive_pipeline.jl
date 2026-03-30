# =============================================================================
#  interactive_pipeline.jl — Step-by-step GRN inference using the public API
#
#  What:
#    Runs the full 6-step mLASSO-StARS pipeline interactively. Each step calls
#    one high-level API function. Intermediate structs remain in the REPL between
#    steps so you can inspect, plot, or debug before continuing.
#
#  Required inputs:
#    geneExprFile        — gene expression matrix (.txt tab-delimited or .arrow)
#    targFile            — target gene list (.txt, one gene per line)
#    regFile             — potential regulator (TF) list (.txt, one TF per line)
#    priorFile           — prior network matrix (.tsv, sparse TF × gene)
#    priorFilePenalties  — prior(s) used to set LASSO penalties (same or different)
#
#  Expected outputs (written under outputDir/<networkBaseName>/):
#    TFA/edges.tsv           — GRN inferred with TFA predictors
#    TFmRNA/edges.tsv        — GRN inferred with TF mRNA predictors
#    Combined/combined_*.tsv — consensus network (max/mean/min aggregation)
#    Combined/TFA/           — refined TFA network re-estimated from consensus prior
#
#  Usage:
#    Step-by-step in a Julia REPL  (recommended for interactive analysis)
#    julia examples/interactive_pipeline.jl
#
#  Compare with:
#    interactive_pipeline_dev.jl  — same steps, module-qualified internal calls
#    run_pipeline.jl              — same pipeline wrapped in a callable function
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================

using InferelatorJL

# =============================================================================
# Configuration — edit these paths and parameters for your dataset
# =============================================================================

outputDir            = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/test"
totSS                = 80
bstarsTotSS          = 5
subsampleFrac        = 0.68
minLambda            = 0.01
maxLambda            = 0.5
totLambdasBstars     = 20
totLambdas           = 40
targetInstability    = 0.05
meanEdgesPerGene     = 20
correlationWeight    = 1
minTargets           = 3
edgeSS               = 0
lambdaBias           = [0.5]
instabilityLevel     = "Network"   # "Network" or "Gene"
useMeanEdgesPerGeneMode = true
combineOpt           = "max"       # "max", "mean", or "min"
zScoreTFA            = true        # z-score targets before TFA estimation
zScoreLASSO          = true        # z-score targets before LASSO regression

geneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_downsample_0.25M.txt"
targFile     = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt"
regFile      = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt"

priorFile          = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"
priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"]
tfaGeneFile        = ""   # optional: restrict TFA estimation to a gene subset

# --- Build output directory name (encodes key run parameters)
subsamplePct    = subsampleFrac * 100
subsampleStr    = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
lambdaStr       = join(replace.(string.(lambdaBias), "." => "p"), "_")
networkBaseName = lowercase(instabilityLevel) * "Lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                  string(meanEdgesPerGene) * "tfsPerGene_" * "subsamplePCT" * subsampleStr
dirOut = joinpath(outputDir, networkBaseName)
mkpath(dirOut)

@info "Configuration" outputDir=dirOut geneExprFile priorFile lambdaBias subsampleFrac

# =============================================================================
# STEP 1 — Load and filter expression data
# =============================================================================
# Loads expression matrix, filters to target genes and potential regulators.
# Inspect: fieldnames(GeneExpressionData), size(data.expressionMat)

data = loadData(geneExprFile, targFile, regFile;
                tfaGeneFile = tfaGeneFile,
                epsilon     = 0.01)

# =============================================================================
# STEP 2 + 3 — Merge degenerate TFs, process prior, estimate TFA
# =============================================================================
# Merges TFs with identical binding profiles, builds the prior matrix,
# and estimates TF activity (TFA) via least-squares.
# Inspect: size(priorData.medTfas), priorData.tfNames

priorData, mergedTFs = loadPrior(data, priorFile; minTargets = minTargets)

estimateTFA(priorData, data;
            edgeSS    = edgeSS,
            zScoreTFA = zScoreTFA,
            outputDir = dirOut)

# =============================================================================
# STEP 4 — Build GRN for each predictor mode
# =============================================================================
# Runs mLASSO-StARS for TFA predictors and TF mRNA predictors separately.
# Outputs instability curves and ranked edge lists to TFA/ and TFmRNA/.

for (tfaMode, modeLabel) in [(true, "TFA"), (false, "TFmRNA")]
    instabilitiesDir = joinpath(dirOut, modeLabel)
    mkpath(instabilitiesDir)

    @info "Building network" mode=modeLabel

    buildNetwork(data, priorData;
                 tfaMode                 = tfaMode,
                 priorFilePenalties      = priorFilePenalties,
                 lambdaBias              = lambdaBias,
                 totSS                   = totSS,
                 bstarsTotSS             = bstarsTotSS,
                 subsampleFrac           = subsampleFrac,
                 minLambda               = minLambda,
                 maxLambda               = maxLambda,
                 totLambdasBstars        = totLambdasBstars,
                 totLambdas              = totLambdas,
                 targetInstability       = targetInstability,
                 meanEdgesPerGene        = meanEdgesPerGene,
                 correlationWeight       = correlationWeight,
                 instabilityLevel        = instabilityLevel,
                 useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode,
                 zScoreLASSO             = zScoreLASSO,
                 outputDir               = instabilitiesDir)
end

# =============================================================================
# STEP 5 — Aggregate TFA + mRNA networks into a consensus network
# =============================================================================
# Combines the two edge lists using max/mean/min stability per (TF, gene) pair.
# Outputs combined_<method>.tsv and combined_<method>_sp.tsv to Combined/.

combinedNetDir = joinpath(dirOut, "Combined")
aggregateNetworks(
    [joinpath(dirOut, "TFA",    "edges.tsv"),
     joinpath(dirOut, "TFmRNA", "edges.tsv")];
    method              = Symbol(combineOpt),
    meanEdgesPerGene    = meanEdgesPerGene,
    useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
    outputDir           = combinedNetDir)

# =============================================================================
# STEP 6 — Re-estimate TFA using the consensus network as a refined prior
# =============================================================================
# Uses the combined network as a new prior to re-estimate TF activity, then
# re-runs mLASSO-StARS. Outputs go to Combined/TFA/.

netsCombinedSparse = joinpath(combinedNetDir, "combined_" * combineOpt * "_sp.tsv")
refineTFA(netsCombinedSparse, data, mergedTFs;
          tfaGeneFile = tfaGeneFile,
          edgeSS      = edgeSS,
          minTargets  = minTargets,
          zScoreTFA   = zScoreTFA,
          exprFile    = geneExprFile,
          targFile    = targFile,
          regFile     = regFile,
          outputDir   = combinedNetDir)

@info "Pipeline complete" outputDir=dirOut
