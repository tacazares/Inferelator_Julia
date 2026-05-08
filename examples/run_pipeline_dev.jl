# =============================================================================
#  run_pipeline_dev.jl — Batch GRN inference via module-qualified internal calls
#
#  What:
#    Identical pipeline to run_pipeline.jl but calls internal functions directly
#    using the InferelatorJL. module prefix. Use this when you need finer control
#    over individual pipeline steps while still running in batch/function mode.
#    All internal functions remain accessible via module-qualified calls even
#    though they are not exported from the package.
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
#    julia examples/run_pipeline_dev.jl
#    or call runInferelator() from another script after: using InferelatorJL
#
#  Compare with:
#    interactive_pipeline_dev.jl  — same internal calls run interactively
#    run_pipeline.jl              — same function using public API
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
    tfaOptions::Vector{String}      = ["", "TFmRNA"],   # "" → TFA, "TFmRNA" → mRNA
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
    modelSelection::String          = "bStARS",   # "bStARS", "EBIC", or "bEBIC"
    gamma::Float64                  = 0.5,         # EBIC sparsity penalty (EBIC / bEBIC only)
    instabilityLevel::String        = "Network",  # "Network" or "Gene"  (bStARS only)
    useMeanEdgesPerGeneMode::Bool   = true,
    combineOpt::String              = "max",       # "max", "mean", or "min"
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

    methodStr = if modelSelection == "bStARS"
        lowercase(instabilityLevel) * "Lambda_" * string(totSS) * "totSS_subsamplePCT" * subsampleStr
    elseif modelSelection == "bEBIC"
        "bEBIC_gamma" * gammaStr * "_" * string(totSS) * "totSS_subsamplePCT" * subsampleStr
    elseif modelSelection == "EBIC"
        "EBIC_gamma" * gammaStr
    end

    networkBaseName = methodStr * "_" * lambdaStr * "_" *
                      string(meanEdgesPerGene) * "tfsPerGene_" *
                      combineOpt * zTFAStr * zLASSOStr * suffix
    dirOut = joinpath(outputDir, networkBaseName)
    mkpath(dirOut)

    # Save all run parameters for reproducibility
    _j(::Nothing)         = "null"
    _j(x::Bool)           = x ? "true" : "false"
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
  "tfaOptions":              $( _j(tfaOptions)              ),
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
    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, geneExprFile)
    InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
    InferelatorJL.loadPotentialRegulators!(data, regFile)
    InferelatorJL.processTFAGenes!(data, tfaGeneFile; outputDir = dirOut)

    # Step 2 — Merge degenerate TFs
    mergedTFsData = mergedTFsResult()
    InferelatorJL.mergeDegenerateTFs(mergedTFsData, priorFile; fileFormat = 2)

    # Step 3 — Process prior and estimate TFA
    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                     mergedTFsData = mergedTFsData,
                                     minTargets    = minTargets)
    InferelatorJL.calculateTFA!(tfaData, data;
                                 edgeSS    = edgeSS,
                                 zTarget   = zScoreTFA,
                                 outputDir = dirOut)

    # Step 3b — Time-lag correction (no-op when timeLagFile = "")
    InferelatorJL.applyTimeLag!(tfaData, data, timeLagFile, timeLag)

    # Step 4 — Build GRN for each predictor mode
    for tfaOpt in tfaOptions
        instabilitiesDir = tfaOpt == "" ? joinpath(dirOut, "TFA") : joinpath(dirOut, "TFmRNA")
        mkpath(instabilitiesDir)

        @info "Building network" tfaOpt=(isempty(tfaOpt) ? "TFA" : tfaOpt)

        grnData  = GrnData()
        buildGrn = BuildGrn()
        InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = tfaOpt)
        InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                             priorFilePenalties = priorFilePenalties,
                                             lambdaBias         = lambdaBias,
                                             tfaOpt             = tfaOpt)

        if modelSelection == "bStARS"
            InferelatorJL.constructSubsamples(data, grnData; totSS = bstarsTotSS, subsampleFrac = subsampleFrac)
            InferelatorJL.bstarsWarmStart(data, tfaData, grnData;
                                           minLambda         = minLambda,
                                           maxLambda         = maxLambda,
                                           totLambdasBstars  = totLambdasBstars,
                                           targetInstability = targetInstability,
                                           zTarget           = zScoreLASSO)
            InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac)
            InferelatorJL.bstartsEstimateInstability(grnData;
                                                      totLambdas        = totLambdas,
                                                      instabilityLevel  = instabilityLevel,
                                                      zTarget           = zScoreLASSO,
                                                      targetInstability = targetInstability,
                                                      outputDir         = instabilitiesDir)
            # Plot on main thread — PyCall crashes if called from a worker thread
            InferelatorJL.plotInstabilityCurves(grnData;
                                                 mode              = :network,
                                                 targetInstability = targetInstability,
                                                 outputDir         = instabilitiesDir)
            GC.gc()   # force Python object cleanup on main thread before next GC-triggering allocation
            InferelatorJL.chooseLambda!(grnData, buildGrn;
                                         instabilityLevel  = instabilityLevel,
                                         targetInstability = targetInstability)

        elseif modelSelection == "EBIC"
            InferelatorJL.ebicSelect!(grnData, buildGrn; gamma = gamma, zScoreLASSO = zScoreLASSO)

        elseif modelSelection == "bEBIC"
            # Per-gene EBIC-optimal lambda on each subsample → lambdaSS, networkStability.
            # bebicFinalFit! is optional (post-hoc only) and not called here.
            InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac)
            InferelatorJL.bebicEstimateLambdas!(grnData, buildGrn;
                                                gamma       = gamma,
                                                minLambda   = minLambda,
                                                maxLambda   = maxLambda,
                                                totLambdas  = totLambdas,
                                                zScoreLASSO = zScoreLASSO,
                                                outputDir   = instabilitiesDir)
        end

        InferelatorJL.rankEdges!(data, tfaData, grnData, buildGrn;
                                    mergedTFsData           = mergedTFsData,
                                    useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode,
                                    meanEdgesPerGene        = meanEdgesPerGene,
                                    correlationWeight       = correlationWeight,
                                    outputDir               = instabilitiesDir)
        writeNetworkTable!(buildGrn; outputDir = instabilitiesDir)
    end

    # Step 5 — Aggregate TFA + mRNA networks into a consensus network
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

    # Step 6 — Re-estimate TFA using the consensus network as a refined prior
    netsCombinedMatrix = joinpath(combinedNetDir, "combined_" * combineOpt * ".tsv")
    InferelatorJL.refineTFA(data, mergedTFsData;
                            priorFile    = netsCombinedMatrix,
                            tfaGeneFile  = tfaGeneFile,
                            edgeSS       = edgeSS,
                            minTargets   = minTargets,
                            zTarget      = zScoreTFA,
                            geneExprFile = geneExprFile,
                            targFile     = targFile,
                            regFile      = regFile,
                            outputDir    = combinedNetDir)

    @info "Pipeline complete" outputDir=dirOut
end

# =============================================================================
# Run — replace paths with your own data
# =============================================================================
# define paths as variables so they can be reused in Step 7 below

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
# ADDED: step 7 for post-pipeline network evaluation using internal computePR directly
using OrderedCollections
import InferelatorJL: computePR

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

        res = computePR(gsFile, outNetFile;
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
