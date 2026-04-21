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
    zScoreLASSO::Bool               = true,
    timeLagFile::String             = "",    # path to 4-column time-lag TSV; leave "" to skip
    timeLag::Real                   = 0.0    # lag value in same units as timeLagFile (e.g. hours)
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

    # Step 3b — Time-lag correction (no-op when timeLagFile = "")
    applyTimeLag(priorData, data, timeLagFile, timeLag)  # ADDED: step 3b, wired from function signature

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
                     mergedTFsData           = mergedTFs,
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
    netsCombinedMatrix = joinpath(combinedNetDir, "combined_" * combineOpt * ".tsv")
    refineTFA(netsCombinedMatrix, data, mergedTFs;
              tfaGeneFile = tfaGeneFile,
              edgeSS      = edgeSS,
              minTargets  = minTargets,
              zScoreTFA   = zScoreTFA,
              timeLagFile = timeLagFile,   # ADDED: forward time-lag to step 6
              timeLag     = timeLag,       # ADDED: forward time-lag to step 6
              exprFile    = geneExprFile,
              targFile    = targFile,
              regFile     = regFile,
              outputDir   = combinedNetDir)

    @info "Pipeline complete" outputDir=dirOut
end

# =============================================================================
# Run — replace paths with your own data
# =============================================================================
# CHANGED: define paths as variables so they can be reused in Step 7 below

geneExprFile       = "/path/to/expression.txt"         # genes × samples (.txt or .arrow)
targFile           = "/path/to/target_genes.txt"       # one gene per line
regFile            = "/path/to/potential_regs.txt"     # one TF per line
priorFile          = "/path/to/prior.tsv"
priorFilePenalties = ["/path/to/prior.tsv"]
outputDir          = "/path/to/output"
combineOpt         = "max"   # must match combineOpt inside runInferelator (default "max")

runInferelator(;
    geneExprFile       = geneExprFile,
    targFile           = targFile,
    regFile            = regFile,
    priorFile          = priorFile,
    priorFilePenalties = priorFilePenalties,
    outputDir          = outputDir
)

# =============================================================================
# STEP 7 — Evaluate networks against a gold standard (optional)
# =============================================================================
# Run after the pipeline completes. Edit gsParam to point to your gold-standard file(s).
# dirOut below reproduces the name built inside runInferelator() with default parameters
# (subsampleFrac=0.68, lambdaBias=[0.5], totSS=80, meanEdgesPerGene=20, instabilityLevel="Network").
# Adjust the networkBaseName string if you changed any of those defaults.
# See examples/plotPR.jl to generate PR curve plots from the saved results.
# ADDED: step 7 for post-pipeline network evaluation
using OrderedCollections

dirOut         = joinpath(outputDir, "networkLambda0p5_80totSS_20tfsPerGene_subsamplePCT68")
combinedNetDir = joinpath(dirOut, "Combined")

outNetFiles = OrderedDict(
    "TFA"      => joinpath(dirOut, "TFA",    "edges_subset.tsv"),
    "TFmRNA"   => joinpath(dirOut, "TFmRNA", "edges_subset.tsv"),
    "Combined" => joinpath(combinedNetDir, "combined_" * combineOpt * "_sp.tsv")
)

# Gold standard(s): name => file path
gsParam = OrderedDict(
    "gsName" => "/path/to/goldStandard.tsv",   # replace with your gold-standard file
)

prFilesByGS = OrderedDict{String, OrderedDict{String, Any}}()
for (legendLabel, outNetFile) in outNetFiles
    @info "Processing network" network=legendLabel file=outNetFile
    filepath = dirname(outNetFile)

    for (gsName, gsFile) in gsParam
        dirPR = joinpath(filepath, "PR", gsName)
        mkpath(dirPR)
        @info "Using GS" gs=gsName saveDir=dirPR

        res = evaluateNetwork(outNetFile, gsFile;
                              gsRegsFile   = regFile,
                              targGeneFile = targFile,
                              breakTies    = true,
                              doPerTF      = false,
                              saveDir      = dirPR)

        if !haskey(prFilesByGS, gsName)
            prFilesByGS[gsName] = OrderedDict{String, Any}()
        end
        prFilesByGS[gsName][legendLabel] = haskey(res, :savedFile) ? res[:savedFile] : res
    end
end

@info "Evaluation complete — see PR/ subdirectories for saved metrics"
