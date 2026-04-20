# =============================================================================
#  Types.jl — All struct definitions for InferelatorJL
#
#  Included FIRST in InferelatorJL.jl so every function file can freely
#  reference any type without worrying about include order.
# =============================================================================


# ── Data loading ──────────────────────────────────────────────────────────────

mutable struct GeneExpressionData
    cellLabels::Vector{String}
    geneNames::Vector{String}
    geneExpressionMat::Matrix{Float64}
    potRegMatmRNA::Matrix{Float64}
    potRegs::Vector{String}
    potRegsmRNA::Vector{String}
    targGenes::Vector{String}
    targGeneMat::Matrix{Float64}
    tfaGenes::Vector{String}
    tfaGeneMat::Matrix{Float64}

    function GeneExpressionData()
        return new(
            [],
            [],
            Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0),
            [],
            [],
            [],
            Matrix{Float64}(undef, 0, 0),
            [],
            Matrix{Float64}(undef, 0, 0)
        )
    end
end


# ── Prior / TF merging ────────────────────────────────────────────────────────

mutable struct mergedTFsResult
    mergedPrior::Union{DataFrame, Nothing}
    mergedTFMap::Union{Matrix{String}, Nothing}

    function mergedTFsResult()
        return new(nothing, nothing)
    end
end


# ── TFA ───────────────────────────────────────────────────────────────────────

mutable struct PriorTFAData
    pRegs::Vector{String}
    pTargs::Vector{String}
    priorMatrix::Matrix{Float64}
    pRegsNoTfa::Vector{String}
    pTargsNoTfa::Vector{String}
    priorMatrixNoTfa::Matrix{Float64}
    noPriorRegs::Vector{String}
    noPriorRegsMat::Matrix{Float64}
    targExpression::Matrix{Float64}
    medTfas::Matrix{Float64}

    function PriorTFAData()
        return new(
            [],
            [],
            Matrix{Float64}(undef, 0, 0),
            [],
            [],
            Matrix{Float64}(undef, 0, 0),
            [],
            Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0)
        )
    end
end


# ── GRN ───────────────────────────────────────────────────────────────────────

mutable struct GrnData
    predictorMat::Matrix{Float64}
    penaltyMat::Matrix{Float64}
    allPredictors::Vector{String}
    subsamps::Matrix{Int64}
    responseMat::Matrix{Float64}
    maxLambdaNet::Float64
    minLambdaNet::Float64
    minLambdas::Matrix{Float64}
    maxLambdas::Matrix{Float64}
    netInstabilitiesUb::Vector{Float64}
    netInstabilitiesLb::Vector{Float64}
    instabilitiesUb::Matrix{Float64}
    instabilitiesLb::Matrix{Float64}
    netInstabilities::Vector{Float64}
    geneInstabilities::Matrix{Float64}
    lambdaRange::Vector{Float64}
    lambdaRangeGene::Vector{Vector{Float64}}
    lambdaRangeWarm::Vector{Float64}             # ADDED: λ grid from warm-start pass (needed for instability bound plots)
    stabilityMat::Array{Float64}
    priorMatProcessed::Matrix{Float64}
    betas::Array{Float64,3}
    targGenes::Vector{String}                    # ADDED: target gene names — row index of stabilityMat (needed for per-gene plots)
    trainInds::Vector{Int}                       # ADDED: training sample indices (excludes leave-out set); needed for R² prediction

    function GrnData()
        return new(
            Matrix{Float64}(undef, 0, 0),        # predictorMat
            Matrix{Float64}(undef, 0, 0),        # penaltyMat
            [],                                   # allPredictors
            Matrix{Int64}(undef, 0, 0),          # subsamps (correctly Int64)
            Matrix{Float64}(undef, 0, 0),        # responseMat
            0.0,                                  # maxLambdaNet
            0.0,                                  # minLambdaNet
            Matrix{Float64}(undef, 0, 0),        # minLambdas
            Matrix{Float64}(undef, 0, 0),        # maxLambdas
            [],                                   # netInstabilitiesUb
            [],                                   # netInstabilitiesLb
            Matrix{Float64}(undef, 0, 0),        # instabilitiesUb
            Matrix{Float64}(undef, 0, 0),        # instabilitiesLb
            [],                                   # netInstabilities
            Matrix{Float64}(undef, 0, 0),        # geneInstabilities
            [],                                   # lambdaRange
            Vector{Vector{Float64}}(undef, 0),   # lambdaRangeGene
            [],                                   # lambdaRangeWarm
            Array{Float64}(undef, 0, 0, 0),      # stabilityMat (3-D)
            Matrix{Float64}(undef, 0, 0),        # priorMatProcessed
            Array{Float64,3}(undef, 0, 0, 0),   # betas
            [],                                   # targGenes
            Int[]                                 # trainInds
        )
    end
end


mutable struct BuildGrn
    networkStability::Matrix{Float64}
    lambda::Union{Float64, Vector{Float64}}
    targs::Vector{String}
    regs::Vector{String}
    rankings::Vector{Float64}
    signedQuantile::Vector{Float64}
    partialCorrelation::Vector{Float64}
    inPrior::Vector{String}
    networkMat::Matrix{Any}
    networkMatSubset::Matrix{Any}
    inPriorVec::Vector{Float64}
    betas::Matrix{Float64}
    mergeTfLocVec::Vector{Float64}

    function BuildGrn()
        return new(
            Matrix{Float64}(undef, 0, 0),  # networkStability
            0.0,                           # lambda
            [],                            # targs
            [],                            # regs
            [],                            # rankings
            [],                            # signedQuantile
            [],                            # partialCorrelation
            [],                            # inPrior
            Matrix{Float64}(undef, 0, 0),  # networkMat
            Matrix{Float64}(undef, 0, 0),  # networkMatSubset
            [],                            # inPriorVec
            Matrix{Float64}(undef, 0, 0),  # betas
            Float64[]                      # mergeTfLocVec
        )
    end
end
