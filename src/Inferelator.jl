module InferelatorJL

# =========================
# Include source files
# =========================
include("00_Data/geneExpression.jl")
include("Utils/dataUtils.jl")
include("01_Prior/mergeDegenerateTFs.jl")
include("00_Data/priorTFA.jl")
include("Utils/networkIO.jl")
include("02_GRN/GRN.jl")

# include("Utils/partialCorrelation.jl")

# =========================
# Helper function to include all files in a folder
# =========================
# function include_all(dir::String)
#     for file in sort(filter(f -> endswith(f, ".jl"), readdir(dir, join=true)))
#         include(file)
#     end
# end

# # ----------------------------------
# # Include submodules
# # ----------------------------------

# # Base data and utils
# include_all(joinpath(@__DIR__, "00_Data"))
# include_all(joinpath(@__DIR__, "Utils"))
# include_all(joinpath(@__DIR__, "01_Prior"))
# include_all(joinpath(@__DIR__, "02_GRN"))

# ---------------------------------------------------------------------------------
# Bring submodules into module namespace
# ---------------------------------------------------------------------------------
using .Data
using .DataUtils
using .MergeDegenerate
using .PriorTFA
using .NetworkIO
using .GRN

# ---------------------------
# Public API
# ---------------------------
export 
    # Core structs
    GeneExpressionData,
    mergedTFsResult,
    PriorTFAData,
    GrnData,
    BuildGrn,

    # Data loading
    loadExpressionData!,
    loadAndFilterTargetGenes!,
    loadPotentialRegulators!,
    processTFAGenes!,

    # Prior / TFA
    # MergeDegenerate,
    processPriorFile!,
    mergeDegenerateTFs,
    calculateTFA!,

    # Data utilities
    convertToLong,
    convertToWide,
    frobeniusNormalize,
    completeDF,
    mergeDFs,
    check_column_norms,
    writeTSVWithEmptyFirstHeader,
    binarizeNumeric!,

    # GRN building
    preparePredictorMat!,
    preparePenaltyMatrix!,
    constructSubsamples,
    bstarsWarmStart,
    bstartsEstimateInstability,
    chooseLambda!,
    rankEdges!,

    # Output
    writeNetworkTable!,
    combineGRNs,
    combineGRNS2,

    #I/O
    saveData, writeNetworkTable!

end