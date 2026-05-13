module InferelatorJL

# ── Dependencies ──────────────────────────────────────────────────────────────
using ArgParse, Arrow, CSV, CategoricalArrays, Colors
using DataFrames, Distributions, FileIO, GLMNet
using InlineStrings, Interpolations, JLD2, Measures
using NamedArrays, OrderedCollections, ProgressMeter
ENV["MPLBACKEND"] = "Agg"   # force non-interactive backend before PyPlot loads
                            # prevents GUI window objects from entering Julia's GC queue,
                            # which would otherwise crash under Threads.@threads (no GIL on worker threads)
using PyPlot, Random, SparseArrays, SpecialFunctions, StatsBase, TickTock

# ── Structs (always first — all type definitions in one place) ────────────────
include("Types.jl")                      # GeneExpressionData, mergedTFsResult, PriorTFAData, GrnData, BuildGrn

# ── Functions ─────────────────────────────────────────────────────────────────
include("data/GeneExpressionData.jl")    # data loading functions
include("data/PseudobulkData.jl")        # pseudobulk aggregation, normalization, batch correction
include("prior/MergeDegenerateTFs.jl")   # TF merging functions
include("data/PriorTFAData.jl")          # prior processing + TFA functions

# ── Core pipeline ─────────────────────────────────────────────────────────────
include("utils/DataUtils.jl")            # Data transformation utilities
include("utils/NetworkIO.jl")            # Save/write network files
include("utils/PartialCorrelation.jl")   # Partial correlation via precision matrix

include("grn/PrepareGRN.jl")            # preparePredictorMat!, preparePenaltyMatrix!, constructSubsamples
include("grn/BuildGRN.jl")              # bstarsWarmStart, bstartsEstimateInstability, chooseLambda!, rankEdges!
include("grn/EBICSelect.jl")            # computeEBIC, ebicSelect!, bebicEstimateLambdas!, bebicFinalFit!
include("grn/AggregateNetworks.jl")     # combineGRNs / aggregateNetworks
include("grn/RefineTFA.jl")             # combineGRNS2 / refineTFA
include("grn/UtilsGRN.jl")             # GRN utility helpers
include("grn/PlotInstability.jl")       # ADDED: plotInstabilityCurves — λ selection diagnostic plots

# ── Metrics ───────────────────────────────────────────────────────────────────
include("metrics/Constants.jl")
include("metrics/MetricUtils.jl")
include("metrics/CalcPR.jl")
include("metrics/CalcR2Pred.jl")
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

    # Pseudobulk input preparation (pure Julia — no R dependency)
    pseudoBulk,
    normalizeMatrix,
    removeBatchEffect,
    aggregateReplicates,
    buildDesignMatrix,

    # Data utilities
    matrixToEdgeList,
    edgeListToMatrix,
    frobeniusNormalize,
    completeDF,
    mergeDFs,
    check_column_norms,
    writeTSVWithEmptyFirstHeader,
    binarizeNumeric!,

    # I/O
    saveData,
    writeNetworkTable!,
    writeBEBICLambdaTable!,
    writeEBICLambdaTable!,
    # R² prediction
    calcR2predFromStabilities,

    # PR metrics
    extractMetricsAtLimit,
    saveSummaryTables,
    computeAUPR,
    computeAUROC

    # -------------------------------------------------------------------------
    # Internal pipeline functions are intentionally NOT exported.
    # They remain accessible via InferelatorJL.<function> or
    # import InferelatorJL: <function> if needed.
    #   loadExpressionData!   loadAndFilterTargetGenes!   loadPotentialRegulators!
    #   processTFAGenes!      processPriorFile!            mergeDegenerateTFs
    #   calculateTFA!         applyTimeLag!                preparePredictorMat!
    #   preparePenaltyMatrix! constructSubsamples          bstarsWarmStart
    #   bstartsEstimateInstability  chooseLambda!          rankEdges!
    #   ebicSelect!           bebicEstimateLambdas!        bebicFinalFit!
    #   computeEBIC
    #   computePR             plotPRCurves                 plotAUPR
    #   loadPRData            plotInstabilityCurves
    # -------------------------------------------------------------------------

end # module InferelatorJL
