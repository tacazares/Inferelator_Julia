module GRN

    include("utilsGRN.jl")

    using ..Data
    using ..DataUtils
    using ..PriorTFA
    using ..MergeDegenerate
    using ..NetworkIO

    using GLMNet
    using Random, Statistics, StatsBase
    using Distributions
    using DataFrames, CategoricalArrays, SparseArrays, CSV, DelimitedFiles
    using LinearAlgebra
    using TickTock
    using JLD2
    using ProgressBars
    using Printf, Dates
    using PyPlot
    using Statistics
    using CSV
    using NamedArrays
    using ArgParse

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
        stabilityMat::Array{Float64}
        priorMatProcessed::Matrix{Float64}
        betas::Array{Float64,3} 
        function GrnData()
            return new(
                Matrix{Float64}(undef, 0, 0), # predictorMat
                Matrix{Float64}(undef, 0, 0), # penaltyMat
                [],                           # allPredictors
                Matrix{Int64}(undef, 0, 0),   # subsamps
                Matrix{Int64}(undef, 0, 0),   # responseMat
                0.0,                           # maxLambdasNet
                0.0,                           # minLambdasNet
                Matrix{Int64}(undef, 0, 0),   # minLambdas
                Matrix{Int64}(undef, 0, 0),   # maxLambdas
                [],                           # netInstabilitiesUb
                [],                           # netInstabilitiesLb
                Matrix{Int64}(undef, 0, 0),   # instabilitiesUb
                Matrix{Int64}(undef, 0, 0),   # instabilitiesLb
                [],                           # netInstabilities
                Matrix{Int64}(undef, 0, 0),   # geneInstabilities
                [],                           # lambdaRange
                Vector{Vector{Float64}}(undef, 0),   # lambdaRangesGene
                Matrix{Int64}(undef, 0, 0),    # stabilityMat
                Matrix{Float64}(undef, 0, 0),  # priorMatProcessed
                Array{Float64,3}(undef, 0, 0, 0) # betas
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
        function BuildGrn()
            return new(
                Matrix{Float64}(undef, 0, 0),   # networkStability
                0.0,                            # lambda
                [],                             # targs
                [],                             # regs
                [],                             # rankings
                [],                             # signedQuantile
                [],                             # partialCorrelation 
                [],                             # inPrior
                Matrix{Float64}(undef, 0, 0),   # networkMat
                Matrix{Float64}(undef, 0, 0),   # networkMatSubset
                [],                             # inPriorVec
                Matrix{Float64}(undef, 0, 0),    # betas
                []                              # mergeTfLocVec
            )
        end
    end

    # Export only the functions that users need
    export preparePredictorMat!, preparePenaltyMatrix!, constructSubsamples, bstarsWarmStart, bstartsEstimateInstability, 
           BuildGrn, GrnData, chooseLambda!, rankEdges!, combineGRNs, combineGRNS2

    include("prepareGRN.jl")   # functions that prepare predictor matrices etc.
    include("buildGRN.jl")     # functions that use GrnData / BuildGrn
    include("aggregateNetworks.jl")   # combineGRNs
    include("refineTFA.jl")  # combineGRNS2


end