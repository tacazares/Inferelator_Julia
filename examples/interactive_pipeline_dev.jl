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
using OrderedCollections
using Dates
# import InferelatorJL: computePR

# =============================================================================
# Configuration — edit these paths and parameters for your dataset
# =============================================================================

outputDir            = "/path/to/output"
tfaOptions           = ["", "TFmRNA"]   # "" → TFA mode, "TFmRNA" → mRNA mode
combineOpt           = "max"            # consensus aggregation method: "max", "mean", or "min" stability per edge
suffix               = ""               # optional label appended to output dir, e.g. "_5KBTSS", "_SCENICprior"

# --- Model selection
modelSelection       = "bStARS"    # "bStARS" : stability-based (default, recommended)
                                   # "EBIC"   : Extended BIC — fast, no subsampling required
                                   # "bEBIC"  : bootstrap EBIC — subsampled, produces selection-frequency scores
gamma                = 0.5         # EBIC sparsity hyperparameter: 0 = BIC, 1 = maximum penalty
                                   # ignored when modelSelection = "bStARS"

# --- Subsampling and stability (bStARS / bEBIC)
totSS                = 80          # total bootstrap subsamples for stability estimation
bstarsTotSS          = 5           # subsamples for warm-start lambda range search (coarser, faster)
subsampleFrac        = 0.68        # fraction of samples per subsample (0.63–0.68 typical)
targetInstability    = 0.05        # instability threshold for bStARS lambda selection (0.05 = 5%)
instabilityLevel     = "Network"   # "Network": one shared lambda for all genes
                                   # "Gene"   : per-gene lambda — slower, more flexible

# --- Lambda grid
minLambda            = 0.01        # LASSO lambda search range — lower bound
maxLambda            = 0.5         # LASSO lambda search range — upper bound
totLambdasBstars     = 20          # lambdas tested in warm-start phase
totLambdas           = 40          # lambdas tested in full stability estimation

# --- Network structure
meanEdgesPerGene     = 20          # average TF regulators retained per target gene
useMeanEdgesPerGeneMode = true     # true : keep meanEdgesPerGene × nGenes edges total
                                   # false: keep all edges above the instability threshold
correlationWeight    = 1           # weight of partial correlation in edge ranking; 0 = stability score only
lambdaBias           = [0.5]       # prior penalty weight(s): 0 = ignore prior, 1 = full prior penalty
                                   # pass multiple values e.g. [0.25, 0.5, 1.0] to sweep

# --- Data processing
minTargets           = 3           # minimum targets a TF must regulate in the prior to be retained
edgeSS               = 0           # TFA edge subsampling replicates; 0 = no subsampling
zScoreTFA            = true        # z-score target expression before TFA estimation
zScoreLASSO          = true        # z-score target expression before LASSO regression

geneExprFile = "/path/to/expression.txt"         # genes × samples (.txt or .arrow)
targFile     = "/path/to/target_genes.txt"       # one gene per line
regFile      = "/path/to/potential_regs.txt"     # one TF per line

priorFile          = "/path/to/prior.tsv"
priorFilePenalties = ["/path/to/prior.tsv"]
tfaGeneFile        = ""   # optional: restrict TFA estimation to a gene subset

# --- Build output directory name (encodes key run parameters)
subsamplePct    = subsampleFrac * 100
subsampleStr    = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
lambdaStr       = join(replace.(string.(lambdaBias), "." => "p"), "_")
zTFAStr         = zScoreTFA   ? "" : "_noZscoreTFA"
zLASSOStr       = zScoreLASSO ? "" : "_noZscoreLASSO"
networkBaseName = lowercase(instabilityLevel) * "Lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                  string(meanEdgesPerGene) * "tfsPerGene_" * "subsamplePCT" * subsampleStr *
                  "_" * combineOpt * zTFAStr * zLASSOStr * suffix
dirOut = joinpath(outputDir, networkBaseName)
mkpath(dirOut)

# Save all run parameters for reproducibility
open(joinpath(dirOut, "run_params.json"), "w") do io
    write(io, """{
  "timestamp":               "$(Dates.now())",
  "geneExprFile":            "$(geneExprFile)",
  "targFile":                "$(targFile)",
  "regFile":                 "$(regFile)",
  "priorFile":               "$(priorFile)",
  "lambdaBias":              $(lambdaBias),
  "totSS":                   $(totSS),
  "bstarsTotSS":             $(bstarsTotSS),
  "subsampleFrac":           $(subsampleFrac),
  "targetInstability":       $(targetInstability),
  "meanEdgesPerGene":        $(meanEdgesPerGene),
  "instabilityLevel":        "$(instabilityLevel)",
  "useMeanEdgesPerGeneMode": $(useMeanEdgesPerGeneMode),
  "combineOpt":              "$(combineOpt)",
  "zScoreTFA":               $(zScoreTFA),
  "zScoreLASSO":             $(zScoreLASSO)
}""")
end

@info "Configuration" outputDir=dirOut geneExprFile=geneExprFile priorFile=priorFile lambdaBias=lambdaBias subsampleFrac=subsampleFrac

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
# STEP 3b — Apply time-lag correction (OPTIONAL — skip if not a time-series)
# =============================================================================
# Adjusts TFA and TF mRNA to account for the delay between TF gene expression
# and active protein. Requires a 4-column TSV (SampleQ, TimeQ, SampleP, TimeP).
# Leave timeLagFile = "" to skip this step entirely.

timeLagFile = "path/to/timeLag.tsv"   # set to "" to skip
timeLag     = 0.5                      # lag in same time units as timeLagFile
InferelatorJL.applyTimeLag!(tfaData, data, timeLagFile, timeLag)

# =============================================================================
# STEP 4 — Build GRN for each predictor mode
# =============================================================================
# Runs lambda selection and edge ranking for each predictor mode.
# The lambda selection branch is controlled by modelSelection in config.

for tfaOpt in tfaOptions
    instabilitiesDir = tfaOpt == "" ? joinpath(dirOut, "TFA") : joinpath(dirOut, "TFmRNA")
    mkpath(instabilitiesDir)

    @info "Building network" mode=(isempty(tfaOpt) ? "TFA" : tfaOpt) modelSelection=modelSelection

    grnData = GrnData()
    InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = tfaOpt)
    InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                         priorFilePenalties = priorFilePenalties,
                                         lambdaBias         = lambdaBias,
                                         tfaOpt             = tfaOpt)
    buildGrn = BuildGrn()

    if modelSelection == "bStARS"
        # Warm-start pass: coarse lambda range on a small number of subsamples
        InferelatorJL.constructSubsamples(data, grnData; totSS = bstarsTotSS, subsampleFrac = subsampleFrac)
        InferelatorJL.bstarsWarmStart(data, tfaData, grnData;
                                       minLambda         = minLambda,
                                       maxLambda         = maxLambda,
                                       totLambdasBstars  = totLambdasBstars,
                                       targetInstability = targetInstability,
                                       zTarget           = zScoreLASSO)
        # Full estimation pass: fine lambda grid on full subsample set
        InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac)
        InferelatorJL.bstartsEstimateInstability(grnData;
                                                  totLambdas       = totLambdas,
                                                  instabilityLevel = instabilityLevel,
                                                  zTarget          = zScoreLASSO,
                                                  outputDir        = instabilitiesDir)
        InferelatorJL.chooseLambda!(grnData, buildGrn;
                                     instabilityLevel  = instabilityLevel,
                                     targetInstability = targetInstability)

    elseif modelSelection == "EBIC"
        # Single full-data LASSO fit per gene — no subsampling required
        InferelatorJL.ebicSelect!(grnData, buildGrn;
                                   gamma       = gamma,
                                   zScoreLASSO = zScoreLASSO)

    elseif modelSelection == "bEBIC"
        # Subsampled EBIC — produces selection-frequency scores like bStARS
        InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac)
        InferelatorJL.bebicSelect!(grnData, buildGrn;
                                    gamma       = gamma,
                                    zScoreLASSO = zScoreLASSO)

    else
        error("modelSelection must be \"bStARS\", \"EBIC\", or \"bEBIC\". Got: \"$modelSelection\"")
    end

    InferelatorJL.rankEdges!(data, tfaData, grnData, buildGrn;
                              mergedTFsData           = mergedTFsData,
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

# =============================================================================
# STEP 7 — Evaluate networks against a gold standard
# =============================================================================
# Compare inferred networks against a gold-standard interaction set.
# Computes precision-recall (PR) metrics and saves results to PR/ subdirectories.
# Edit gsParam to point to your gold-standard file(s).
# See examples/plotPR.jl to generate PR curve plots from the saved results.

outNetFiles = OrderedDict(
    "TFA"      => joinpath(dirOut, "TFA",    "edges_subset.tsv"),
    "TFmRNA"   => joinpath(dirOut, "TFmRNA", "edges_subset.tsv"),
    "Combined" => joinpath(combinedNetDir, "combined_" * combineOpt * "_sp.tsv")
)

# Gold standard(s): name => file path
gsParam = OrderedDict(
    "gsName" => "/path/to/goldStandard.tsv",   # replace with your gold-standard name and file
)

prFilesByGS = OrderedDict{String, OrderedDict{String, Any}}()
for (legendLabel, outNetFile) in outNetFiles
    @info "Processing network" network=legendLabel file=outNetFile
    filepath = dirname(outNetFile)

    for (gsName, gsFile) in gsParam
        dirPR = joinpath(filepath, "PR", gsName)
        mkpath(dirPR)
        @info "Using GS" gs=gsName saveDir=dirPR

        res = InferelatorJL.computePR(gsFile, outNetFile;
                        gsRegsFile   = regFile,
                        targGeneFile = targFile,
                        breakTies    = true,
                        doPerTF      = false,
                        xLimitRecall = 0.1,
                        saveDir      = dirPR)

        if !haskey(prFilesByGS, gsName)
            prFilesByGS[gsName] = OrderedDict{String, Any}()
        end
        prFilesByGS[gsName][legendLabel] = haskey(res, :savedFile) ? res[:savedFile] : res
    end
end

@info "Evaluation complete — see PR/ subdirectories for saved metrics"


