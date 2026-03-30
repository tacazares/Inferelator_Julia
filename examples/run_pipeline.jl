# =============================================================================
#  run_pipeline.jl — Batch GRN inference wrapped in a callable function (public API)
#
#  What:
#    Wraps the full 6-step mLASSO-StARS pipeline in runInferelator(), which
#    accepts all parameters as keyword arguments with sensible defaults.
#    Use this for batch execution, scripted cluster jobs, or when running
#    multiple parameter sweeps programmatically.
#    All steps call high-level public API functions (no internal functions exposed).
#
#  Required inputs:
#    geneExprFile        — gene expression matrix (.txt tab-delimited or .arrow)
#    targFile            — target gene list (.txt, one gene per line)
#    regFile             — potential regulator (TF) list (.txt, one TF per line)
#    priorFile           — prior network matrix (.tsv, sparse TF × gene)
#    priorFilePenalties  — prior(s) used to set LASSO penalties (same or different)
#    outputDir           — root directory for all outputs
#
#  Expected outputs (written under outputDir/<networkBaseName>/):
#    TFA/edges.tsv           — GRN inferred with TFA predictors
#    TFmRNA/edges.tsv        — GRN inferred with TF mRNA predictors
#    Combined/combined_*.tsv — consensus network (max/mean/min aggregation)
#    Combined/TFA/           — refined TFA network re-estimated from consensus prior
#
#  Usage:
#    julia examples/run_pipeline.jl
#    or call runInferelator() from another script after: using InferelatorJL
#
#  Compare with:
#    interactive_pipeline.jl  — same steps run interactively (not wrapped)
#    run_pipeline_dev.jl      — same function using module-qualified internal calls
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================

using InferelatorJL

# =============================================================================
# Pipeline function
# =============================================================================

function runInferelator(;
    geneExprFile::String,
    targFile::String,
    regFile::String,
    priorFile::String,
    priorFilePenalties::Vector{String},
    tfaGeneFile::String             = "",
    outputDir::String,
    totSS::Int                      = 80,
    bstarsTotSS::Int                = 5,
    subsampleFrac::Float64          = 0.68,
    minLambda::Float64              = 0.01,
    maxLambda::Float64              = 0.5,
    totLambdasBstars::Int           = 20,
    totLambdas::Int                 = 40,
    targetInstability::Float64      = 0.05,
    meanEdgesPerGene::Int           = 20,
    correlationWeight::Int          = 1,
    minTargets::Int                 = 3,
    edgeSS::Int                     = 0,
    lambdaBias::Vector{Float64}     = [0.5],
    instabilityLevel::String        = "Network",  # "Network" or "Gene"
    useMeanEdgesPerGeneMode::Bool   = true,
    combineOpt::String              = "max",       # "max", "mean", or "min"
    zScoreTFA::Bool                 = true,
    zScoreLASSO::Bool               = true
)
    # Build output directory name (encodes key run parameters)
    subsamplePct    = subsampleFrac * 100
    subsampleStr    = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
    lambdaStr       = join(replace.(string.(lambdaBias), "." => "p"), "_")
    networkBaseName = lowercase(instabilityLevel) * "Lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                      string(meanEdgesPerGene) * "tfsPerGene_" * "subsamplePCT" * subsampleStr
    dirOut = joinpath(outputDir, networkBaseName)
    mkpath(dirOut)

    @info "Starting pipeline" outputDir=dirOut geneExprFile priorFile lambdaBias subsampleFrac

    # Step 1 — Load and filter expression data
    data = loadData(geneExprFile, targFile, regFile;
                    tfaGeneFile = tfaGeneFile,
                    epsilon     = 0.01)

    # Steps 2 + 3 — Merge degenerate TFs, process prior, estimate TFA
    priorData, mergedTFs = loadPrior(data, priorFile; minTargets = minTargets)

    estimateTFA(priorData, data;
                edgeSS    = edgeSS,
                zScoreTFA = zScoreTFA,
                outputDir = dirOut)

    # Step 4 — Build GRN for each predictor mode
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

    # Step 5 — Aggregate TFA + mRNA networks into a consensus network
    combinedNetDir = joinpath(dirOut, "Combined")
    aggregateNetworks(
        [joinpath(dirOut, "TFA",    "edges.tsv"),
         joinpath(dirOut, "TFmRNA", "edges.tsv")];
        method              = Symbol(combineOpt),
        meanEdgesPerGene    = meanEdgesPerGene,
        useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
        outputDir           = combinedNetDir)

    # Step 6 — Re-estimate TFA using the consensus network as a refined prior
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
end

# =============================================================================
# Run — replace paths with your own data
# =============================================================================

runInferelator(
    geneExprFile       = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt",
    targFile           = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt",
    regFile            = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt",
    priorFile          = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv",
    priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"],
    outputDir          = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/test"
)
