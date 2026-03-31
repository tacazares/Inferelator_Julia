# =============================================================================
#  interactive_pipeline_dev.jl — Step-by-step GRN inference via internal functions
#
#  What:
#    Identical pipeline to interactive_pipeline.jl but calls internal functions
#    directly using the InferelatorJL. module prefix. Use this when you need
#    finer control over individual steps (e.g., inspecting intermediate matrices,
#    swapping in a custom subfunction, or debugging a specific stage).
#    All internal functions remain accessible via module-qualified calls even
#    though they are not exported from the package.
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
#    Step-by-step in a Julia REPL  (recommended for inspecting intermediate state)
#    julia examples/interactive_pipeline_dev.jl
#
#  Compare with:
#    interactive_pipeline.jl  — same steps via high-level public API
#    run_pipeline_dev.jl      — same internal calls wrapped in a callable function
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================
using Revise
using InferelatorJL

# =============================================================================
# Configuration — edit these paths and parameters for your dataset
# =============================================================================

outputDir            = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/test"
tfaOptions           = ["", "TFmRNA"]   # "" → TFA mode, "TFmRNA" → mRNA mode
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
# Inspect: fieldnames(GeneExpressionData), size(data.expressionMat), data.geneNames

data = GeneExpressionData()
InferelatorJL.loadExpressionData!(data, geneExprFile)
InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
InferelatorJL.loadPotentialRegulators!(data, regFile)
InferelatorJL.processTFAGenes!(data, tfaGeneFile; outputDir = dirOut)

# =============================================================================
# STEP 2 — Merge degenerate TFs
# =============================================================================
# Identifies TFs with identical binding profiles in the prior and merges them.
# Inspect: mergedTFsData.mergedTFs, mergedTFsData.tfNames

mergedTFsData = mergedTFsResult()
InferelatorJL.mergeDegenerateTFs(mergedTFsData, priorFile; fileFormat = 2)

# =============================================================================
# STEP 3 — Process prior and estimate TFA
# =============================================================================
# Builds the filtered prior matrix and estimates TF activity via least-squares.
# Inspect: tfaData.priorMat, tfaData.medTfas, size(tfaData.medTfas)

tfaData = PriorTFAData()
InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                 mergedTFsData = mergedTFsData,
                                 minTargets    = minTargets)
InferelatorJL.calculateTFA!(tfaData, data;
                             edgeSS    = edgeSS,
                             zTarget   = zScoreTFA,
                             outputDir = dirOut)

# =============================================================================
# STEP 4 — Build GRN for each predictor mode
# =============================================================================
# Runs subsampling, warm-start lambda selection, instability estimation,
# lambda selection, and edge ranking for each predictor mode.

for tfaOpt in tfaOptions
    instabilitiesDir = tfaOpt == "" ? joinpath(dirOut, "TFA") : joinpath(dirOut, "TFmRNA")
    mkpath(instabilitiesDir)

    @info "Building network" tfaOpt=(isempty(tfaOpt) ? "TFA" : tfaOpt)

    grnData = GrnData()
    InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = tfaOpt)
    InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                         priorFilePenalties = priorFilePenalties,
                                         lambdaBias         = lambdaBias,
                                         tfaOpt             = tfaOpt)
    InferelatorJL.constructSubsamples(data, grnData; totSS = bstarsTotSS, subsampleFrac = subsampleFrac)
    InferelatorJL.bstarsWarmStart(data, tfaData, grnData;
                                   minLambda         = minLambda,
                                   maxLambda         = maxLambda,
                                   totLambdasBstars  = totLambdasBstars,
                                   targetInstability = targetInstability,
                                   zTarget           = zScoreLASSO)
    InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac)
    InferelatorJL.bstartsEstimateInstability(grnData;
                                              totLambdas       = totLambdas,
                                              instabilityLevel = instabilityLevel,
                                              zTarget          = zScoreLASSO,
                                              outputDir        = instabilitiesDir)

    buildGrn = BuildGrn()
    InferelatorJL.chooseLambda!(grnData, buildGrn;
                                 instabilityLevel  = instabilityLevel,
                                 targetInstability = targetInstability)
    InferelatorJL.rankEdges!(data, tfaData, grnData, buildGrn;
                              useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode,
                              meanEdgesPerGene        = meanEdgesPerGene,
                              correlationWeight       = correlationWeight,
                              outputDir               = instabilitiesDir)
    writeNetworkTable!(buildGrn; outputDir = instabilitiesDir)
end

# =============================================================================
# STEP 5 — Aggregate TFA + mRNA networks into a consensus network
# =============================================================================
# Combines edge lists using max/mean/min stability per (TF, gene) pair.

combinedNetDir = joinpath(dirOut, "Combined")
nets2combine   = [
    joinpath(dirOut, "TFA",    "edges.tsv"),
    joinpath(dirOut, "TFmRNA", "edges.tsv")
]
InferelatorJL.aggregateNetworks(nets2combine;
                                method              = Symbol(combineOpt),
                                meanEdgesPerGene    = meanEdgesPerGene,
                                useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
                                outputDir           = combinedNetDir)

# =============================================================================
# STEP 6 — Re-estimate TFA using the consensus network as a refined prior
# =============================================================================
# Uses combined network as new prior, re-estimates TFA, re-runs mLASSO-StARS.

netsCombinedSparse = joinpath(combinedNetDir, "combined_" * combineOpt * "_sp.tsv")
InferelatorJL.refineTFA(data, mergedTFsData;
                        priorFile    = netsCombinedSparse,
                        tfaGeneFile  = tfaGeneFile,
                        edgeSS       = edgeSS,
                        minTargets   = minTargets,
                        zTarget      = zScoreTFA,
                        geneExprFile = geneExprFile,
                        targFile     = targFile,
                        regFile      = regFile,
                        outputDir    = combinedNetDir)

@info "Pipeline complete" outputDir=dirOut
