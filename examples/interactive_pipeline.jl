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
using Revise
using InferelatorJL
using Dates

# =============================================================================
# Configuration — edit these paths and parameters for your dataset
# =============================================================================

outputDir            = "/path/to/output"
combineOpt = :max       # consensus aggregation method: "max", "mean", or "min" stability per edge
suffix               = ""          # optional label appended to output dir, e.g. "_5KBTSS", "_SCENICprior"

# --- Model selection ---------------------------------------------------------
modelSelection = :bStARS   # "bStARS" : stability-based (default, recommended)
                             # "EBIC"   : Extended BIC — fast, no subsampling required
                             # "bEBIC"  : bootstrap EBIC — subsampled, selection-frequency scores
gamma          = 0.5         # EBIC sparsity penalty: 0 = BIC, 1 = maximum sparsity
                             # (EBIC / bEBIC only — ignored by bStARS)

# --- Subsampling (bStARS and bEBIC) ------------------------------------------
totSS         = 80    # total bootstrap subsamples
subsampleFrac = 0.68  # fraction of samples per subsample (0.63–0.68 typical)

# --- bStARS-only settings ----------------------------------------------------
bstarsTotSS       = 5     # subsamples for warm-start lambda range search (coarser pass)
targetInstability = 0.05  # instability threshold for lambda selection (5% default)
instabilityLevel = :network  # "Network": one shared lambda across all genes
                               # "Gene"   : per-gene lambda — slower, more flexible

# --- Lambda grid -------------------------------------------------------------
minLambda        = nothing  # lower bound (nothing = auto-computed per gene as max|X'y|/n × eps)
maxLambda        = nothing  # upper bound (nothing = auto-computed per gene as max|X'y|/n)
totLambdasBstars = 20       # grid points in warm-start pass  (bStARS only)
totLambdas       = 100      # grid points in full estimation pass (log-spaced; all methods)

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
zScoreTFA            = false       # z-score target expression before TFA estimation
zScoreLASSO          = false       # z-score target expression before LASSO regression
timeLagFile          = ""          # path to 4-column time-lag TSV; "" skips step 3b
timeLag              = 0.0         # lag value in same units as timeLagFile (e.g. hours)

geneExprFile = "/path/to/expression.txt"         # genes × samples (.txt or .arrow)
targFile     = "/path/to/target_genes.txt"       # one gene per line
regFile      = "/path/to/potential_regs.txt"     # one TF per line

priorFile          = "/path/to/prior.tsv"
priorFilePenalties = ["/path/to/prior.tsv"]
tfaGeneFile        = ""   # optional: restrict TFA estimation to a gene subset

# --- Build output directory name (encodes method and key run parameters)
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

@info "Configuration" outputDir=dirOut geneExprFile=geneExprFile priorFile=priorFile lambdaBias=lambdaBias subsampleFrac=subsampleFrac

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
# STEP 3b — Apply time-lag correction (OPTIONAL — skip if not a time-series)
# =============================================================================
# Adjusts TFA and TF mRNA to account for the delay between TF gene expression
# and active protein. Requires a 4-column TSV (SampleQ, TimeQ, SampleP, TimeP).
# Leave timeLagFile = "" to skip this step entirely.
# ADDED: step 3b for optional time-lag correction

applyTimeLag(priorData, data, timeLagFile, timeLag)

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
                 mergedTFsData           = mergedTFs,
                 modelSelection          = modelSelection,
                 gamma                   = gamma,
                 outputDir               = instabilitiesDir)
end

# =============================================================================
# STEP 5 — Aggregate TFA + mRNA networks into a consensus network
# =============================================================================
# Combines the two edge lists using max/mean/min stability per (TF, gene) pair.
# Outputs combined_<method>_sp.tsv (long/edge list) and combined_<method>.tsv (wide/matrix) to Combined/.

combinedNetDir = joinpath(dirOut, "Combined")
aggregateNetworks(
    [joinpath(dirOut, "TFA",    "edges.tsv"),
     joinpath(dirOut, "TFmRNA", "edges.tsv")];
    method = combineOpt,
    meanEdgesPerGene    = meanEdgesPerGene,
    useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
    outputDir           = combinedNetDir)

# =============================================================================
# STEP 6 — Re-estimate TFA using the consensus network as a refined prior
# =============================================================================
# Uses the combined network as a new prior to re-estimate TF activity, then
# re-runs mLASSO-StARS. Outputs go to Combined/TFA/.

netsCombinedMatrix = joinpath(combinedNetDir, "combined_" * string(combineOpt) * ".tsv")
refineTFA(netsCombinedMatrix, data, mergedTFs;
          tfaGeneFile = tfaGeneFile,
          edgeSS      = edgeSS,
          minTargets  = minTargets,
          zScoreTFA   = zScoreTFA,
          exprFile    = geneExprFile,
          targFile    = targFile,
          regFile     = regFile,
          outputDir   = combinedNetDir)

@info "Pipeline complete" outputDir=dirOut

# =============================================================================
# STEP 7 — Evaluate networks against a gold standard
# =============================================================================
# Compare inferred networks against a gold-standard interaction set.
# Computes precision-recall (PR) metrics and saves results to PR/ subdirectories.
# Edit gsParam to point to your gold-standard file(s).
# See examples/plotPR.jl to generate PR curve plots from the saved results.
# ADDED: step 7 for network evaluation using the public evaluateNetwork API
using OrderedCollections

outNetFiles = OrderedDict(
    "TFA"      => joinpath(dirOut, "TFA",    "edges_subset.tsv"),
    "TFmRNA"   => joinpath(dirOut, "TFmRNA", "edges_subset.tsv"),
    "Combined" => joinpath(combinedNetDir, "combined_" * string(combineOpt) * "_sp.tsv")
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
