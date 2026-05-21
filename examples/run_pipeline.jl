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
using Dates

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
    minLambda::Union{Float64, Nothing} = nothing,
    maxLambda::Union{Float64, Nothing} = nothing,
    totLambdasBstars::Int           = 20,
    totLambdas::Int                 = 100,
    targetInstability::Float64      = 0.05,
    meanEdgesPerGene::Int           = 20,
    correlationWeight::Int          = 1,
    minTargets::Int                 = 3,
    edgeSS::Int                     = 0,
    lambdaBias::Vector{Float64}     = [0.5],
    modelSelection::Symbol          = :bStARS,    # :bStARS, :EBIC, or :bEBIC
    gamma::Float64                  = 0.5,         # EBIC sparsity penalty (EBIC / bEBIC only)
    instabilityLevel::Symbol        = :network,   # :network (single λ) or :gene (per-gene λ)  (bStARS only)
    useMeanEdgesPerGeneMode::Bool   = true,
    combineOpt::Symbol              = :max,        # :max, :mean, or :min
    zScoreTFA::Bool                 = false,
    zScoreLASSO::Bool               = false,
    timeLagFile::String             = "",    # path to 4-column time-lag TSV; leave "" to skip
    timeLag::Real                   = 0.0,   # lag value in same units as timeLagFile (e.g. hours)
    suffix::String                  = ""     # optional custom label appended to output dir, e.g. "_5KBTSS"
)
    # Build output directory name (encodes method and key run parameters)
    subsamplePct = subsampleFrac * 100
    subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
    lambdaStr    = join(replace.(string.(lambdaBias), "." => "p"), "_")
    zTFAStr      = zScoreTFA   ? "" : "_noZscoreTFA"
    zLASSOStr    = zScoreLASSO ? "" : "_noZscoreLASSO"
    gammaStr     = replace(string(gamma), "." => "p")

    methodStr = if modelSelection == :bStARS
        string(instabilityLevel) * "Lambda_" * string(totSS) * "totSS_subsamplePCT" * subsampleStr
    elseif modelSelection == :bEBIC
        "bEBIC_gamma" * gammaStr * "_" * string(totSS) * "totSS_subsamplePCT" * subsampleStr
    elseif modelSelection == :EBIC
        "EBIC_gamma" * gammaStr
    end

    networkBaseName = methodStr * "_" * lambdaStr * "_" *
                      string(meanEdgesPerGene) * "tfsPerGene_" *
                      string(combineOpt) * zTFAStr * zLASSOStr * suffix
    dirOut = joinpath(outputDir, networkBaseName)
    mkpath(dirOut)

    # Save all run parameters for reproducibility
    _j(::Nothing)         = "null"
    _j(x::Bool)           = x ? "true" : "false"
    _j(x::Symbol)         = "\"" * string(x) * "\""
    _j(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
    _j(x::AbstractVector) = "[" * join(_j.(x), ", ") * "]"
    _j(x)                 = string(x)

    open(joinpath(dirOut, "run_params.json"), "w") do io
        write(io, """{
  "timestamp":               "$( Dates.now()               )",
  "geneExprFile":            $( _j(geneExprFile)            ),
  "targFile":                $( _j(targFile)                ),
  "regFile":                 $( _j(regFile)                 ),
  "priorFile":               $( _j(priorFile)               ),
  "priorFilePenalties":      $( _j(priorFilePenalties)      ),
  "tfaGeneFile":             $( _j(tfaGeneFile)             ),
  "modelSelection":          $( _j(modelSelection)          ),
  "gamma":                   $( _j(gamma)                   ),
  "lambdaBias":              $( _j(lambdaBias)              ),
  "minLambda":               $( _j(minLambda)               ),
  "maxLambda":               $( _j(maxLambda)               ),
  "totSS":                   $( _j(totSS)                   ),
  "bstarsTotSS":             $( _j(bstarsTotSS)             ),
  "totLambdasBstars":        $( _j(totLambdasBstars)        ),
  "totLambdas":              $( _j(totLambdas)              ),
  "subsampleFrac":           $( _j(subsampleFrac)           ),
  "targetInstability":       $( _j(targetInstability)       ),
  "meanEdgesPerGene":        $( _j(meanEdgesPerGene)        ),
  "correlationWeight":       $( _j(correlationWeight)       ),
  "instabilityLevel":        $( _j(instabilityLevel)        ),
  "useMeanEdgesPerGeneMode": $( _j(useMeanEdgesPerGeneMode) ),
  "combineOpt":              $( _j(combineOpt)              ),
  "minTargets":              $( _j(minTargets)              ),
  "edgeSS":                  $( _j(edgeSS)                  ),
  "zScoreTFA":               $( _j(zScoreTFA)               ),
  "zScoreLASSO":             $( _j(zScoreLASSO)             ),
  "timeLagFile":             $( _j(timeLagFile)             ),
  "timeLag":                 $( _j(timeLag)                 ),
  "suffix":                  $( _j(suffix)                  )
}""")
    end

    @info "Starting pipeline" outputDir=dirOut geneExprFile=geneExprFile priorFile=priorFile lambdaBias=lambdaBias subsampleFrac=subsampleFrac

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
                     modelSelection          = modelSelection,
                     gamma                   = gamma,
                     outputDir               = instabilitiesDir)
    end

    # Step 5 — Aggregate TFA + mRNA networks into a consensus network
    combinedNetDir = joinpath(dirOut, "Combined")
    aggregateNetworks(
        [joinpath(dirOut, "TFA",    "edges.tsv"),
         joinpath(dirOut, "TFmRNA", "edges.tsv")];
        method = combineOpt,
        meanEdgesPerGene    = meanEdgesPerGene,
        useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
        outputDir           = combinedNetDir)

    # Step 6 — Re-estimate TFA using the consensus network as a refined prior
    netsCombinedMatrix = joinpath(combinedNetDir, "combined_" * string(combineOpt) * ".tsv")
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
combineOpt = :max   # must match combineOpt inside runInferelator (default "max")

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
# (subsampleFrac=0.68, lambdaBias=[0.5], totSS=80, meanEdgesPerGene=20, instabilityLevel = :network).
# Adjust the networkBaseName string if you changed any of those defaults.
# See examples/plotPR.jl to generate PR curve plots from the saved results.
# ADDED: step 7 for post-pipeline network evaluation
using OrderedCollections

dirOut         = joinpath(outputDir, "networkLambda0p5_80totSS_20tfsPerGene_subsamplePCT68")
combinedNetDir = joinpath(dirOut, "Combined")

outNetFiles = OrderedDict(
    "TFA"      => joinpath(dirOut, "TFA",    "edges_subset.tsv"),
    "TFmRNA"   => joinpath(dirOut, "TFmRNA", "edges_subset.tsv"),
    "Combined" => joinpath(combinedNetDir, "combined_" * string(combineOpt) * "_sp.tsv")
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
