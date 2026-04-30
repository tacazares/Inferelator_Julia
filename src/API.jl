# =============================================================================
#  InferelatorJL — Public API
#  src/API.jl
#
#  High-level wrapper functions that compose the internal !-mutating calls
#  into clean, single-responsibility entry points.
#  All heavy lifting stays in the submodule files; this layer only
#  orchestrates and exposes a stable interface.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 1 · Load & filter all expression data
# -----------------------------------------------------------------------------
"""
    loadData(exprFile, targFile, regFile; tfaGeneFile="", epsilon=0.01)

Load and filter all expression inputs into a `GeneExpressionData` struct.

# Arguments
- `exprFile`    : Path to gene expression matrix (TSV or Arrow, genes × samples)
- `targFile`    : Path to target gene list (genes to model as responses)
- `regFile`     : Path to potential regulator list (candidate TFs)
- `tfaGeneFile` : Optional path to gene list used for TFA estimation
- `epsilon`     : Minimum variance threshold for target gene filtering (default 0.01)

# Returns
`GeneExpressionData` with fields populated:
`geneExpressionMat`, `targGeneMat`, `potRegMatmRNA`, `tfaGeneMat`
"""
function loadData(
    exprFile::String,
    targFile::String,
    regFile::String;
    tfaGeneFile::String = "",
    epsilon::Float64    = 0.01
)::GeneExpressionData

    data = GeneExpressionData()
    loadExpressionData!(data, exprFile)
    loadAndFilterTargetGenes!(data, targFile; epsilon = epsilon)
    loadPotentialRegulators!(data, regFile)
    processTFAGenes!(data, tfaGeneFile)
    return data
end


# -----------------------------------------------------------------------------
# STEP 2 + 3 · Merge degenerate TFs then process prior & estimate TFA
# -----------------------------------------------------------------------------
"""
    loadPrior(data, priorFile; minTargets=3, mergeDegenerate=true, connector="_")

Merge degenerate TFs, process the prior matrix, and return both the
`PriorTFAData` struct and the `mergedTFsResult` map.

# Arguments
- `data`            : Populated `GeneExpressionData` from `loadData`
- `priorFile`       : Path to TF × Gene prior matrix (TSV)
- `minTargets`      : Minimum number of targets a TF must have to be retained (default 3)
- `mergeDegenerate` : Whether to collapse TFs with identical target sets (default true)
- `connector`       : String used to join meta-TF names (default "_")

# Returns
Tuple `(PriorTFAData, mergedTFsResult)`
"""
function loadPrior(
    data::GeneExpressionData,
    priorFile::String;
    minTargets::Int      = 3,
    mergeDegenerate::Bool = true,
    connector::String    = "_"
)::Tuple{PriorTFAData, mergedTFsResult}

    mergedTFs = mergedTFsResult()
    if mergeDegenerate
        mergeDegenerateTFs(mergedTFs, priorFile; fileFormat = 2, connector = connector)
    end

    priorData = PriorTFAData()
    processPriorFile!(priorData, data, priorFile;
                      mergedTFsData = mergedTFs,
                      minTargets    = minTargets)
    return priorData, mergedTFs
end


# -----------------------------------------------------------------------------
# STEP 3 (cont.) · Estimate TF activity
# -----------------------------------------------------------------------------
"""
    estimateTFA(priorData, data; edgeSS=0, zScoreTFA=true, outputDir=".")

Estimate TF activity (TFA) by solving `prior * TFA ≈ targetExpression`
via least squares, with optional bootstrap subsampling of targets.

# Arguments
- `priorData`   : `PriorTFAData` from `loadPrior`
- `data`        : `GeneExpressionData` from `loadData`
- `edgeSS`      : Number of edge subsamples (0 = no subsampling)
- `zScoreTFA`   : Z-score target expression before solving TFA (default true)
- `outputDir`   : Directory for intermediate output files (default ".")

# Returns
TFA matrix as `Matrix{Float64}` (TFs × samples), also stored in `priorData.medTfas`
"""
function estimateTFA(
    priorData::PriorTFAData,
    data::GeneExpressionData;
    edgeSS::Int       = 0,
    zScoreTFA::Bool   = true,
    outputDir::String = "."
)::Matrix{Float64}

    calculateTFA!(priorData, data;
                  edgeSS    = edgeSS,
                  zTarget   = zScoreTFA,
                  outputDir = outputDir)
    return priorData.medTfas
end


# -----------------------------------------------------------------------------
# STEP 3b (optional) · Apply time-lag correction to TFA and TF mRNA
# -----------------------------------------------------------------------------
"""
    applyTimeLag(priorData, data, timeLagFile, timeLag)

Apply a time-lag correction to TFA and TF mRNA estimates after `estimateTFA`.

For each consecutive time-point pair (P → Q) in `timeLagFile`, the value at
sample Q is adjusted to account for the expected delay between TF mRNA and
active TF protein:

    adjusted[Q] = original[Q] − (original[Q] − original[P]) / Δt × timeLag

Three matrices are updated in-place:
- `priorData.medTfas`             — prior-based TFA
- `priorData.noPriorRegsMat`      — mRNA for regulators absent from the prior
- `data.potRegMatmRNA`            — mRNA for all potential regulators

# Arguments
- `priorData`    : `PriorTFAData` with `medTfas` already populated by `estimateTFA`
- `data`         : `GeneExpressionData` from `loadData`
- `timeLagFile`  : path to a 4-column tab-separated file (with header).
  Columns: sample Q name, time of Q, sample P name, time of P.
  Only time-series pairs need to be listed; unpaired samples are unchanged.
- `timeLag`      : delay between TF mRNA and protein activity, in the same
  time units as the values in `timeLagFile` (e.g., hours).

# Example
```julia
data, priorData, _ = loadData(...), loadPrior(...), ...
estimateTFA(priorData, data)
applyTimeLag(priorData, data, "timeLag.tsv", 0.5)
```

# Reference
Bonneau et al. (2006) Genome Biology; Miraldi et al. `integratePrior_estTFA_timeLag`.
"""
function applyTimeLag(
    priorData::PriorTFAData,
    data::GeneExpressionData,
    timeLagFile::String,
    timeLag::Real
)
    applyTimeLag!(priorData, data, timeLagFile, timeLag)
end


# -----------------------------------------------------------------------------
# STEP 4 · Build a single GRN (one predictor mode)
# -----------------------------------------------------------------------------
"""
    buildNetwork(data, priorData; kwargs...) → BuildGrn

Run the full mLASSO-StARS pipeline for one predictor mode and return the
ranked edge table.

Internally runs:
`preparePredictorMat!` → `preparePenaltyMatrix!` → `constructSubsamples` →
`bstarsWarmStart` → `constructSubsamples` (full) → `bstartsEstimateInstability` →
`chooseLambda!` → `rankEdges!` → `writeNetworkTable!`

# Arguments
- `data`                  : `GeneExpressionData`
- `priorData`             : `PriorTFAData`
- `tfaMode`               : true = TFA predictors, false = raw mRNA for all TFs
- `priorFilePenalties`    : Prior file(s) used to build the penalty matrix
- `lambdaBias`            : Penalty reduction factor for prior edges (default [0.5])
- `totSS`                 : Total subsamples for fine instability estimation (default 80)
- `bstarsTotSS`           : Subsamples for warm-start (default 5)
- `subsampleFrac`         : Fraction of samples per subsample (default 0.63)
- `minLambda`             : Lower bound of λ search range (default 0.01)
- `maxLambda`             : Upper bound of λ search range (default 0.5)
- `totLambdasBstars`      : λ grid points in warm-start pass (default 20)
- `totLambdas`            : λ grid points in fine estimation pass (default 40)
- `targetInstability`     : Instability threshold for λ selection (default 0.05)
- `meanEdgesPerGene`      : Max edges retained per target gene (default 20)
- `correlationWeight`     : Weight of partial correlation in edge scoring (default 1)
- `instabilityLevel`      : "Network" (single λ) or "Gene" (per-gene λ)
- `useMeanEdgesPerGeneMode`: Enforce per-gene edge cap (default true)
- `zTarget`               : Z-score targets during regression (default true)
- `leaveOutSampleList`    : Path to a text file listing samples to hold out (one per line);
                            held-out samples are excluded from subsampling and used as the
                            test set for R²_pred evaluation. Pass "" to use all samples.
- `mergedTFsData`         : `mergedTFsResult` returned by `loadPrior`; required for TF
                            de-merging after regression. Pass `nothing` to skip expansion.
- `modelSelection`        : Lambda selection method: "bStARS" (default), "EBIC", or "bEBIC".
                            "bStARS" uses stability-based selection (current default).
                            "EBIC" uses a single full-data fit scored by Extended BIC.
                            "bEBIC" uses EBIC on each subsample and takes the median lambda.
- `gamma`                 : EBIC hyperparameter controlling the sparsity penalty (default 0.5).
                            Only used when `modelSelection` is "EBIC" or "bEBIC".
                            gamma=0 reduces to standard BIC; gamma=0.5 is recommended for GRN.
- `outputDir`             : Output directory for edges.tsv and stability arrays

# Returns
`BuildGrn` with `networkStability`, `signedQuantile`, `networkMat`, etc.
"""
function buildNetwork(
    data::GeneExpressionData,
    priorData::PriorTFAData;
    tfaMode::Bool                   = true,
    priorFilePenalties::Vector{String} = String[],
    lambdaBias::Vector{Float64}     = [0.5],
    totSS::Int                      = 80,
    bstarsTotSS::Int                = 5,
    subsampleFrac::Float64          = 0.63,
    minLambda::Float64              = 0.01,
    maxLambda::Float64              = 0.5,
    totLambdasBstars::Int           = 20,
    totLambdas::Int                 = 40,
    targetInstability::Float64      = 0.05,
    meanEdgesPerGene::Int           = 20,
    correlationWeight::Int          = 1,
    instabilityLevel::String        = "Network",
    useMeanEdgesPerGeneMode::Bool   = true,
    zScoreLASSO::Bool               = true,
    leaveOutSampleList::String      = "",
    mergedTFsData::Union{mergedTFsResult, Nothing} = nothing,
    modelSelection::String          = "bStARS",
    gamma::Float64                  = 0.5,
    outputDir::String               = "."
)::BuildGrn

    tfaOpt  = tfaMode ? "" : "TFmRNA"
    loList  = isempty(leaveOutSampleList) ? nothing : leaveOutSampleList

    grnData = GrnData()
    preparePredictorMat!(grnData, data, priorData; tfaOpt = tfaOpt)
    preparePenaltyMatrix!(data, grnData;
                          priorFilePenalties = priorFilePenalties,
                          lambdaBias         = lambdaBias,
                          tfaOpt             = tfaOpt)

    buildGrn = BuildGrn()

    if modelSelection == "bStARS"
        # Warm-start pass (coarse lambda range)
        constructSubsamples(data, grnData;
                            totSS              = bstarsTotSS,
                            subsampleFrac      = subsampleFrac,
                            leaveOutSampleList = loList)
        bstarsWarmStart(data, priorData, grnData;
                        minLambda         = minLambda,
                        maxLambda         = maxLambda,
                        totLambdasBstars  = totLambdasBstars,
                        targetInstability = targetInstability,
                        zTarget           = zScoreLASSO)

        # Fine estimation pass
        constructSubsamples(data, grnData;
                            totSS              = totSS,
                            subsampleFrac      = subsampleFrac,
                            leaveOutSampleList = loList)
        bstartsEstimateInstability(grnData;
                                   totLambdas        = totLambdas,
                                   instabilityLevel  = instabilityLevel,
                                   zTarget           = zScoreLASSO,
                                   targetInstability = targetInstability,
                                   outputDir         = outputDir)
        chooseLambda!(grnData, buildGrn;
                      instabilityLevel  = instabilityLevel,
                      targetInstability = targetInstability)

    elseif modelSelection == "EBIC"
        # Single full-data fit; no subsampling needed
        ebicSelect!(grnData, buildGrn;
                    gamma         = gamma,
                    zScoreLASSO   = zScoreLASSO)

    elseif modelSelection == "bEBIC"
        # Subsampled EBIC — requires subsample indices
        constructSubsamples(data, grnData;
                            totSS              = totSS,
                            subsampleFrac      = subsampleFrac,
                            leaveOutSampleList = loList)
        bebicSelect!(grnData, buildGrn;
                     gamma         = gamma,
                     zScoreLASSO   = zScoreLASSO)

    else
        error("modelSelection must be \"bStARS\", \"EBIC\", or \"bEBIC\". Got: \"$modelSelection\"")
    end

    rankEdges!(data, priorData, grnData, buildGrn;
               mergedTFsData           = mergedTFsData,
               useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode,
               meanEdgesPerGene        = meanEdgesPerGene,
               correlationWeight       = correlationWeight,
               outputDir               = outputDir)
    writeNetworkTable!(buildGrn; outputDir = outputDir)

    return buildGrn
end


# -----------------------------------------------------------------------------
# STEP 6 · Recalculate TFA using the combined network as a refined prior
# -----------------------------------------------------------------------------
"""
    refineTFA(combinedNetFile, data, mergedTFs; kwargs...) → Matrix{Float64}

Use the aggregated consensus network as a new prior and re-estimate TF
activity, yielding a data-driven TFA matrix that reflects the final GRN.

# Arguments
- `combinedNetFile` : Path to sparse combined network TSV (combined_max_sp.tsv)
- `data`            : `GeneExpressionData`
- `mergedTFs`       : `mergedTFsResult` from `loadPrior`
- `tfaGeneFile`  : Optional gene list for TFA (default "")
- `edgeSS`       : Edge subsampling replicates for TFA (default 0)
- `minTargets`   : Minimum targets per TF (default 3)
- `zScoreTFA`    : Z-score target expression before solving TFA (default true)
- `timeLagFile`  : Path to 4-column time-lag TSV; omit or pass "" to skip (default "")
- `timeLag`      : Time-lag value in the same units as `timeLagFile` (default 0.0)
- `exprFile`     : Original expression file path
- `targFile`     : Target gene file path
- `regFile`      : Regulator file path
- `outputDir`    : Directory for refined TFA output (default ".")

# Returns
Refined TFA matrix as `Matrix{Float64}` (regulators × samples)
"""
function refineTFA(
    combinedNetFile::String,
    data::GeneExpressionData,
    mergedTFs::mergedTFsResult;
    tfaGeneFile::String  = "",
    edgeSS::Int          = 0,
    minTargets::Int      = 3,
    zScoreTFA::Bool      = true,
    timeLagFile::String  = "",
    timeLag::Real        = 0.0,
    exprFile::String     = "",
    targFile::String     = "",
    regFile::String      = "",
    outputDir::String    = "."
)

    # Dispatch to internal refineTFA(data::GeneExpressionData, ...) — different first-arg type
    refineTFA(data, mergedTFs;
              priorFile    = combinedNetFile,
              tfaGeneFile  = tfaGeneFile,
              edgeSS       = edgeSS,
              minTargets   = minTargets,
              zTarget      = zScoreTFA,
              timeLagFile  = timeLagFile,
              timeLag      = timeLag,
              geneExprFile = exprFile,
              targFile     = targFile,
              regFile      = regFile,
              outputDir    = outputDir)
end


# -----------------------------------------------------------------------------
# STEP 7 · Evaluate a network against a gold standard
# -----------------------------------------------------------------------------
"""
    evaluateNetwork(networkFile, goldStandard; kwargs...) → Dict

Evaluate a ranked edge list against a gold-standard interaction set using
precision-recall (PR) metrics.  Wraps the internal `computePR` function.

# Arguments
- `networkFile`       : Path to inferred GRN edge list (TSV with TF, Gene, score columns)
- `goldStandard`      : Path to gold-standard network TSV (sparse TF × Gene format)
- `gsRegsFile`        : Optional regulator list to restrict evaluation to shared TFs
                        (default "" → use all TFs present in both files)
- `targGeneFile`      : Optional target gene list to restrict the evaluation universe
                        (default "" → use all genes present in both files)
- `breakTies`         : Randomly break ties in edge scores before ranking (default true)
- `partialAUPRlimit`  : Recall limit for partial AUPR calculation (default 0.1)
- `doPerTF`           : Compute per-TF PR curves in addition to global PR (default false)
- `saveDir`           : Directory to write PR result files; "" skips saving (default "")

# Returns
`Dict` with keys `:aupr`, `:auroc`, `:precisions`, `:recalls`, `:randPR`, and
optionally `:savedFile` (path to the saved results file when `saveDir` is set).

# Example
```julia
res = evaluateNetwork(
    joinpath(dirOut, "TFA", "edges_subset.tsv"),
    "/path/to/goldStandard.tsv";
    gsRegsFile   = regFile,
    targGeneFile = targFile,
    saveDir      = joinpath(dirOut, "TFA", "PR", "ChIP")
)
println("AUPR: ", res[:aupr])
```
"""
function evaluateNetwork(
    networkFile::String,
    goldStandard::String;
    gsRegsFile::String        = "",
    targGeneFile::String      = "",
    breakTies::Bool           = true,
    partialAUPRlimit::Float64 = 0.1,
    doPerTF::Bool             = false,
    saveDir::String           = ""
)
    # CHANGED: replaced incorrect computeMacroMetrics call with computePR
    computePR(goldStandard, networkFile;
              gsRegsFile       = isempty(gsRegsFile)   ? nothing : gsRegsFile,
              targGeneFile     = isempty(targGeneFile)  ? nothing : targGeneFile,
              breakTies        = breakTies,
              partialAUPRlimit = partialAUPRlimit,
              doPerTF          = doPerTF,
              saveDir          = isempty(saveDir)       ? nothing : saveDir)
end


# -----------------------------------------------------------------------------
# CONVENIENCE · Full pipeline in one call
# -----------------------------------------------------------------------------
"""
    inferGRN(exprFile, targFile, regFile, priorFile; outputDir="results", kwargs...)

Run the complete InferelatorJL pipeline end-to-end:

1. Load & filter expression data
2. Merge degenerate TFs + process prior
3. Estimate TFA
4. Build TFA-mode network
5. Build mRNA-mode network
6. Aggregate both networks
7. Recalculate TFA on combined network

All keyword arguments are forwarded to the relevant sub-functions.

# Key keyword arguments
- `outputDir::String = "results"` — root directory for all outputs
- `totSS::Int = 80` — total bootstrap subsamples
- `subsampleFrac::Float64 = 0.63` — fraction of samples per subsample
- `lambdaBias::Vector{Float64} = [0.5]` — penalty reduction for prior edges
- `targetInstability::Float64 = 0.05` — instability threshold for λ selection
- `meanEdgesPerGene::Int = 20` — maximum retained edges per target gene
- `timeLagFile::String = ""` — path to 4-column TSV for time-lag correction; `""` skips step 3b
- `timeLag::Real = 0.0` — lag value in same time units as `timeLagFile`
- `leaveOutSampleList::String = ""` — path to a text file (one sample name per line) listing
  samples to exclude from subsampling; the held-out samples form the test set used by
  `calcR2predFromStabilities` for out-of-sample R² evaluation
- `modelSelection::String = "bStARS"` — lambda selection method: "bStARS", "EBIC", or "bEBIC"
- `gamma::Float64 = 0.5` — EBIC sparsity hyperparameter; only used when modelSelection != "bStARS"

# Returns
`BuildGrn` from the TFA-mode network (combined results written to `outputDir`)
"""
function inferGRN(
    exprFile::String,
    targFile::String,
    regFile::String,
    priorFile::String;
    outputDir::String               = "results",
    tfaGeneFile::String             = "",
    epsilon::Float64                = 0.01,
    minTargets::Int                 = 3,
    edgeSS::Int                     = 0,
    zScoreTFA::Bool                 = true,
    zScoreLASSO::Bool               = true,
    priorFilePenalties::Vector{String} = String[],
    lambdaBias::Vector{Float64}     = [0.5],
    totSS::Int                      = 80,
    bstarsTotSS::Int                = 5,
    subsampleFrac::Float64          = 0.63,
    minLambda::Float64              = 0.01,
    maxLambda::Float64              = 0.5,
    totLambdasBstars::Int           = 20,
    totLambdas::Int                 = 40,
    targetInstability::Float64      = 0.05,
    meanEdgesPerGene::Int           = 20,
    correlationWeight::Int          = 1,
    instabilityLevel::String        = "Network",
    useMeanEdgesPerGeneMode::Bool   = true,
    combineMethod::Symbol           = :max,
    timeLagFile::String             = "",
    timeLag::Real                   = 0.0,
    leaveOutSampleList::String      = "",
    modelSelection::String          = "bStARS",
    gamma::Float64                  = 0.5
)::BuildGrn

    mkpath(outputDir)
    tfaDir   = joinpath(outputDir, "TFA")
    mRNADir  = joinpath(outputDir, "TFmRNA")
    combDir  = joinpath(outputDir, "Combined")
    mkpath(tfaDir);  mkpath(mRNADir);  mkpath(combDir)

    # Shared kwargs for buildNetwork
    netKwargs = (
        priorFilePenalties      = isempty(priorFilePenalties) ? [priorFile] : priorFilePenalties,
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
        modelSelection          = modelSelection,
        gamma                   = gamma,
    )

    # Steps 1–3
    data                = loadData(exprFile, targFile, regFile;
                                   tfaGeneFile = tfaGeneFile, epsilon = epsilon)
    priorData, mergedTFs = loadPrior(data, priorFile; minTargets = minTargets)
    estimateTFA(priorData, data; edgeSS = edgeSS, zScoreTFA = zScoreTFA,
                outputDir = outputDir)
    if !isempty(timeLagFile)
        applyTimeLag(priorData, data, timeLagFile, timeLag)
    end

    # Step 4 — TFA mode
    tfaGrn  = buildNetwork(data, priorData; tfaMode = true,
                            leaveOutSampleList = leaveOutSampleList,
                            mergedTFsData      = mergedTFs,
                            netKwargs..., outputDir = tfaDir)

    # Step 4 — mRNA mode
    mrnaGrn = buildNetwork(data, priorData; tfaMode = false,
                            leaveOutSampleList = leaveOutSampleList,
                            mergedTFsData      = mergedTFs,
                            netKwargs..., outputDir = mRNADir)

    # Step 5 — aggregate
    aggregateNetworks(
        [joinpath(tfaDir,  "edges.tsv"),
         joinpath(mRNADir, "edges.tsv")];
        method              = combineMethod,
        meanEdgesPerGene    = meanEdgesPerGene,
        useMeanEdgesPerGene = useMeanEdgesPerGeneMode,
        outputDir           = combDir
    )

    # Step 6 — refine TFA
    combinedMatrix = joinpath(combDir, "combined_" * string(combineMethod) * ".tsv")
    refineTFA(combinedMatrix, data, mergedTFs;
              tfaGeneFile = tfaGeneFile,
              edgeSS      = edgeSS,
              minTargets  = minTargets,
              zScoreTFA   = zScoreTFA,
              timeLagFile = timeLagFile,
              timeLag     = timeLag,
              exprFile    = exprFile,
              targFile    = targFile,
              regFile     = regFile,
              outputDir   = combDir)

    return tfaGrn
end
