module InferelatorJL

# ── Dependencies ──────────────────────────────────────────────────────────────
using ArgParse, Arrow, CSV, CategoricalArrays, Colors
using DataFrames, Distributions, FileIO, GLMNet
using InlineStrings, Interpolations, JLD2, Measures
using NamedArrays, OrderedCollections, ProgressBars
using PyPlot, Random, SparseArrays, StatsBase, TickTock

# ── Structs (always first — all type definitions in one place) ────────────────
include("Types.jl")                      # GeneExpressionData, mergedTFsResult, PriorTFAData, GrnData, BuildGrn

# ── Functions ─────────────────────────────────────────────────────────────────
include("data/GeneExpressionData.jl")    # data loading functions
include("prior/MergeDegenerateTFs.jl")   # TF merging functions
include("data/PriorTFAData.jl")          # prior processing + TFA functions

# ── Core pipeline ─────────────────────────────────────────────────────────────
include("utils/DataUtils.jl")            # Data transformation utilities
include("utils/NetworkIO.jl")            # Save/write network files
include("utils/PartialCorrelation.jl")   # Partial correlation via precision matrix

include("grn/PrepareGRN.jl")            # preparePredictorMat!, preparePenaltyMatrix!, constructSubsamples
include("grn/BuildGRN.jl")              # bstarsWarmStart, bstartsEstimateInstability, chooseLambda!, rankEdges!
include("grn/AggregateNetworks.jl")     # combineGRNs / aggregateNetworks
include("grn/RefineTFA.jl")             # combineGRNS2 / refineTFA
include("grn/UtilsGRN.jl")             # GRN utility helpers

# ── Metrics ───────────────────────────────────────────────────────────────────
include("metrics/Constants.jl")
include("metrics/MetricUtils.jl")
include("metrics/CalcPR.jl")
include("metrics/Metrics.jl")
include("metrics/plotting/PlotSingle.jl")
include("metrics/plotting/PlotBatch.jl")

# ── Public API ────────────────────────────────────────────────────────────────
include("API.jl")

# ── Exports ───────────────────────────────────────────────────────────────────
export
    # Structs
    GeneExpressionData,
    mergedTFsResult,
    PriorTFAData,
    GrnData,
    BuildGrn,

    # High-level API  (src/API.jl)
    loadData,
    loadPrior,
    estimateTFA,
    applyTimeLag,
    buildNetwork,
    aggregateNetworks,
    refineTFA,
    evaluateNetwork,
    inferGRN,

    # Data utilities
    convertToLong,
    convertToWide,
    frobeniusNormalize,
    completeDF,
    mergeDFs,
    check_column_norms,
    writeTSVWithEmptyFirstHeader,
    binarizeNumeric!,

    # I/O
    saveData,
    writeNetworkTable!

    # -------------------------------------------------------------------------
    # Internal pipeline functions are intentionally NOT exported.
    # They remain accessible via InferelatorJL.<function> or
    # import InferelatorJL: <function> if needed.
    #   loadExpressionData!   loadAndFilterTargetGenes!   loadPotentialRegulators!
    #   processTFAGenes!      processPriorFile!            mergeDegenerateTFs
    #   calculateTFA!         applyTimeLag!                preparePredictorMat!
    #   preparePenaltyMatrix! constructSubsamples          bstarsWarmStart
    #   bstartsEstimateInstability  chooseLambda!          rankEdges!
    #   computePR             plotPRCurves                 plotAUPR
    #   loadPRData
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

end # module InferelatorJL
