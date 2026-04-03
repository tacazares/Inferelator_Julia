"""
GRN

Functions:
- preparePredictorMat!
- preparePenaltyMatrix!
- constructSubsamples
- bstarsWarmStart
- bstartsEstimateInstability

Dependencies:
"""

# mutable struct GrnData
#     predictorMat::Matrix{Float64}
#     penaltyMat::Matrix{Float64}
#     allPredictors::Vector{String}
#     subsamps::Matrix{Int64}
#     responseMat::Matrix{Float64}
#     maxLambdaNet::Float64
#     minLambdaNet::Float64
#     minLambdas::Matrix{Float64}
#     maxLambdas::Matrix{Float64}
#     netInstabilitiesUb::Vector{Float64}
#     netInstabilitiesLb::Vector{Float64}
#     instabilitiesUb::Matrix{Float64}
#     instabilitiesLb::Matrix{Float64}
#     netInstabilities::Vector{Float64}
#     geneInstabilities::Matrix{Float64}
#     lambdaRange::Vector{Float64}
#     lambdaRangeGene::Vector{Vector{Float64}}
#     stabilityMat::Array{Float64}
#     priorMatProcessed::Matrix{Float64}
#     betas::Array{Float64,3} 
#     function GrnData()
#         return new(
#             Matrix{Float64}(undef, 0, 0), # predictorMat
#             Matrix{Float64}(undef, 0, 0), # penaltyMat
#             [],                           # allPredictors
#             Matrix{Int64}(undef, 0, 0),   # subsamps
#             Matrix{Int64}(undef, 0, 0),   # responseMat
#             0.0,                           # maxLambdasNet
#             0.0,                           # minLambdasNet
#             Matrix{Int64}(undef, 0, 0),   # minLambdas
#             Matrix{Int64}(undef, 0, 0),   # maxLambdas
#             [],                           # netInstabilitiesUb
#             [],                           # netInstabilitiesLb
#             Matrix{Int64}(undef, 0, 0),   # instabilitiesUb
#             Matrix{Int64}(undef, 0, 0),   # instabilitiesLb
#             [],                           # netInstabilities
#             Matrix{Int64}(undef, 0, 0),   # geneInstabilities
#             [],                           # lambdaRange
#             Vector{Vector{Float64}}(undef, 0),   # lambdaRangesGene
#             Matrix{Int64}(undef, 0, 0),    # stabilityMat
#             Matrix{Float64}(undef, 0, 0),  # priorMatProcessed
#             Array{Float64,3}(undef, 0, 0, 0) # betas
#         )
#     end
# end

function preparePredictorMat!(grnData::GrnData, expressionData::GeneExpressionData, priorData::PriorTFAData; tfaOpt::String = "")
    if tfaOpt != ""
        @info "TF mRNA mode (no TFA)"
        pRegs = priorData.pRegsNoTfa;
        pTargs = priorData.pTargsNoTfa;
        priorMatrix = priorData.priorMatrixNoTfa;
    else
        pRegs = priorData.pRegs;
        pTargs = priorData.pTargs;
        priorMatrix = priorData.priorMatrix;
    end
    # TFs in potRegs that have expression data that arent in prior. Get the index for these TFs 
    # in the TF expression matrix
    uniNoPriorRegs = setdiff(expressionData.potRegsmRNA, pRegs)
    uniNoPriorRegInds = findall(in(uniNoPriorRegs), expressionData.potRegsmRNA)

    # allPredictors include both the pRegs and the potRegs that wernt in prior
    allPredictors = vcat(pRegs, uniNoPriorRegs)
    totPreds = length(allPredictors)

    # Create new prior matrix that contains target genes in the same order as targGenes. If using
    # TFA, missing TFs not in prior that we have expression data for
    targGeneInds = findall(in(pTargs), expressionData.targGenes)
    priorGeneInds = findall(in(expressionData.targGenes), pTargs)
    totTargGenes = length(expressionData.targGenes)
    totPRegs = length(pRegs)
    priorMat = zeros(totTargGenes,totPreds)
    priorMat[targGeneInds,1:totPRegs] = priorMatrix[priorGeneInds,:]

    # # predictorMat will be TFA (when available) and TFmRNA when TFA not available    
    # predictorMat = [priorData.medTfas; expressionData.potRegMatmRNA[uniNoPriorRegInds,:]]
    # # If not using TFA, just set predictorMat to mRNA
    # if tfaOpt != "" # use the mRNA levels of TFs
    #     currPredMat = zeros(totPreds,length(expressionData.cellLabels))
    #     for prend = 1:totPreds
    #         prendInd = findall(x -> x==allPredictors[prend],expressionData.potRegsmRNA)
    #         currPredMat[prend,:] = expressionData.potRegMatmRNA[prendInd,:]
    #     end            
    #     predictorMat = currPredMat
    #     println("TF mRNA used.")
    # end

    # predictorMat will be TFA (when available) and TFmRNA when TFA not available   
    if tfaOpt == ""
        # TFA: stack TFA estimates on top of mRNA for TFs not in prior
        predictorMat = [priorData.medTfas; expressionData.potRegMatmRNA[uniNoPriorRegInds,:]]
    else
        # If not using TFA: use mRNA for all predictors
        currPredMat = zeros(totPreds, length(expressionData.cellLabels))
        for prend = 1:totPreds
            prendInd = findall(x -> x == allPredictors[prend], expressionData.potRegsmRNA)
            currPredMat[prend, :] = expressionData.potRegMatmRNA[prendInd, :]
        end
        predictorMat = currPredMat
        @info "TF mRNA used as predictors" 
    end
    grnData.predictorMat = predictorMat
    grnData.allPredictors = allPredictors
    grnData.responseMat = expressionData.targGeneMat
    grnData.priorMatProcessed = priorMat
    grnData.targGenes = copy(expressionData.targGenes)   # ADDED: store gene names for per-gene instability plots
end


function preparePenaltyMatrix!(expressionData::GeneExpressionData, grnData::GrnData;
                               priorFilePenalties::Vector{String} = String[],
                               lambdaBias::Vector{Float64}        = [0.5],
                               tfaOpt::String                     = "")
    #1. Update Penalty Matrix
    # Create a dictionary to store the minimum lambda for each interaction
    minLambdaDict = Dict{Tuple{String,String}, Float64}()
    penaltyMatrix = ones(length(expressionData.targGenes),length(grnData.allPredictors))

    # Iterate through each prior file and its associated lambda
    for (filePath, lambda) in zip(priorFilePenalties, lambdaBias)
        # Read the prior file
        priorData = readdlm(filePath)
        # priorData = CSV.read(filePath, DataFrame; delim = "\t")

        # Extract gene names and TF names from the prior file
        priorGenes = priorData[2:end, 1]
        # Get TF names, filtering out any empty or missing entries
        priorTFs = filter(tf -> !ismissing(tf) && !isempty(string(tf)), priorData[1, :])

        # Create indices mapping for faster lookup
        geneToIdx = Dict(gene => i for (i, gene) in enumerate(expressionData.targGenes))
        tfToIdx = Dict(tf => i for (i, tf) in enumerate(grnData.allPredictors))
        # Process the interactions
        for (i, gene) in enumerate(priorGenes)
            for (j, tf) in enumerate(priorTFs)
                # println("  Trying ($gene, $tf): ")
                if priorData[i+1, j+1] != 0 && haskey(geneToIdx, gene) && haskey(tfToIdx, tf)
                    interaction = (gene, tf)
                    if !haskey(minLambdaDict, interaction) || lambda < minLambdaDict[interaction]
                        minLambdaDict[interaction] = lambda
                    end
                end
            end
        end
    end

    # Apply the penalties using the minimum lambda values
    for ((gene, tf), minLambda) in minLambdaDict
        geneIdx = findfirst(==(gene), expressionData.targGenes)
        tfIdx = findfirst(==(tf), grnData.allPredictors)
        penaltyMatrix[geneIdx, tfIdx] = minLambda
    end

    # 2. 
    totPreds = length(grnData.allPredictors)
    if tfaOpt !== ""
        ## set lambda penalty to infinity for positive feedback edges where TF 
        # mRNA levels serves both as gene expression and TFA estimate
        for pr = 1:totPreds   # Changed this from length(expressionData.potRegs) to length(grnData.allPredictors)
            targInd = findall(x -> x==grnData.allPredictors[pr], expressionData.targGenes)
            if length(targInd) > 0 # set lambda penalty to infinity, avoid predicting a TF's mRNA based on its own mRNA level
                penaltyMatrix[targInd,pr] .= Inf # i.e., target gene is its own predictor
            end
        end    
    else # have to set prior inds to zero for TFs in TFA that don't have prior info
        for pr = 1:totPreds    # Changed this from expressionData.potRegs to grnData.allPredictors  
            if sum(abs.(grnData.priorMatProcessed[:,pr])) == 0 # we have no target edges to estimate TF's TFA
                targInd = findall(x -> x==grnData.allPredictors[pr], expressionData.targGenes)
                if length(targInd) > 0 # And TF is in the predictor set
                    penaltyMatrix[targInd,pr] .= Inf
                end
            end
        end
    end

    grnData.penaltyMat = penaltyMatrix
end


function constructSubsamples(expressionData::GeneExpressionData, grnData::GrnData; 
                            leaveOutSampleList::Union{Vector{Vector{String}}, Nothing}=nothing,totSS = 200, subsampleFrac = 0.63)
    totSamps = length(expressionData.cellLabels)

    if !(leaveOutSampleList in (nothing, ""))
        @info "Leave-out set detected" leaveOutSampleList 
        # get leave-out set of samples
        fin = open(leaveOutSampleList)
        C = readdlm(fin,skipstart=0)
        C = convert(Matrix{String}, C)
        close(fin)
        testInds = Tuple.(findall(in(C), expressionData.cellLabels))
        testInds = first.(testInds)
        trainInds = setdiff(1:totSamps,testInds)
    else        
        @info "Full gene expression matrix used (no leave-out set)"
        trainInds = 1:totSamps # all training samples used
        testInds = []
    end

    subsampleSize = floor(subsampleFrac*length(trainInds))
    subsampleSize = convert(Int64, subsampleSize)
    # get subsamp indices
    subsamps = zeros(totSS,subsampleSize)
    for ss = 1:totSS
        randSubs = randperm(totSamps)
        randSubs = randSubs[1:subsampleSize]
        subsamps[ss,:] = randSubs
    end
    subsamps = convert(Matrix{Int}, subsamps)
    grnData.subsamps = subsamps
end

function bstarsWarmStart(expressionData::GeneExpressionData, tfaData::PriorTFAData, grnData::GrnData; minLambda = 0.01, maxLambda = 0.5, totLambdasBstars = 20, totSS = 5, targetInstability = 0.05, zTarget = false)
    # Determine the lambda levels to test
    lambdaRange = collect(range(minLambda, stop = maxLambda, length = totLambdasBstars))
    #lamLog10step = 1/totLambdasBstars
    #logLamRange =  log10(minLambda):lamLog10step:log10(maxLambda)
    #lambdaRange = 10 .^ (logLamRange)
    lambdaRange = reverse(lambdaRange)
    totResponses = size(grnData.responseMat)[1]

    instabilitiesLb = zeros(totResponses,totLambdasBstars)
    instabilitiesUb = zeros(totResponses,totLambdasBstars)
    minLambdas = zeros(totResponses,1)
    maxLambdas = zeros(totResponses,1)

    responsePredInds = Vector{Vector{Int}}(undef,0)
    for res = 1:totResponses
        currWeights = grnData.penaltyMat[res,:]
        push!(responsePredInds,findall(x -> x!=Inf, currWeights))
    end
    
    netInstabilitiesLb = zeros(totResponses, totLambdasBstars)
    netInstabilitiesUb = zeros(totResponses, totLambdasBstars)
    totEdges = zeros(totResponses)   # denominator for network Instabilities

    theta2save = Vector{Matrix{Float64}}(undef, totResponses)  # Debugging #

    Threads.@threads for res in ProgressBar(1:totResponses) # can be a parfor loop
        predInds = responsePredInds[res]
        currPredNum = length(predInds)
        penaltyFactor = grnData.penaltyMat[res, predInds]
        totEdges[res] = currPredNum
        ssVals = zeros(totLambdasBstars,currPredNum)
        for ss = 1:totSS
            subsamp = grnData.subsamps[ss,:]
            dt = fit(ZScoreTransform, grnData.predictorMat[predInds, subsamp], dims=2)
            currPreds = transpose(StatsBase.transform(dt, grnData.predictorMat[predInds, subsamp]))
            if zTarget
                dt = fit(ZScoreTransform, grnData.responseMat[res, subsamp], dims=1)
                currResponses = StatsBase.transform(dt, grnData.responseMat[res, subsamp])
            else
                currResponses = grnData.responseMat[res, subsamp]
            end
            lsoln = glmnet(currPreds, currResponses, penalty_factor = penaltyFactor, lambda = lambdaRange, alpha = 1.0)
            currBetas = lsoln.betas # flip so that the lambdas are increasing
            # ssVals += abs.(sign.(currBetas))'
            ssVals .+= abs.(sign.(currBetas))'
        end
        theta2 = (1/totSS)*ssVals # empirical edge probability
        theta2save[res] = theta2  #For  Debugging #
        instabilitiesLb[res,:] = 2 * (1/currPredNum) .* sum(theta2 .* (1 .- theta2), dims=2) # bStARS lower bound
        netInstabilitiesLb[res,:] = currPredNum*(instabilitiesLb[res,:])
        theta2mean = sum(theta2,dims=2)./currPredNum
        instabilitiesUb[res,:] = 2 * theta2mean .* (1 .- theta2mean) # bStARS upper bound
        netInstabilitiesUb[res,:] = currPredNum*instabilitiesUb[res,:]
    end

    totEdges = sum(totEdges)
    netInstabilitiesLb = sum(netInstabilitiesLb, dims=1)[:]
    netInstabilitiesUb = sum(netInstabilitiesUb, dims=1)[:]

    for res = 1:totResponses
        # take the supremum, find max Lambda, and set all smaller lambdas equal to that value 
        maxLb = findmax(instabilitiesLb[res,:])
        maxLbInd = findall(x -> x == maxLb[1], instabilitiesLb[res,:])
        maxLb = maxLb[1]
        instabilitiesLb[res,maxLbInd[end]:end] .= maxLb
        maxUb = findmax(instabilitiesUb[res,:])
        maxUbInd = findall(x -> x == maxUb[1], instabilitiesUb[res,:])
        maxUb = maxUb[1]
        instabilitiesUb[res,maxUbInd[end]:end] .= maxUb
        # find the minimum lambda for the gene, based on maximum for upper bound
        # we are less interested in high instability lambdas, so okay to use
        # upper bound
        xx = findmin(abs.(instabilitiesLb[res,:] .- targetInstability))
        #xx = findall(x -> x == xx[1],abs.(instabilitiesLb[res,:] .- targetMaxInstability))
        xx = xx[2]
        minLambdas[res] = (lambdaRange)[xx[end]] # to the right
        # find the lambda nearest the min instability worth considering, use
        # upper bound as that will be sure to find an lambda >= target instability lambda 
        xx = findmin(abs.(instabilitiesUb[res,:] .- targetInstability))
        xx = findall(x -> x == xx[1], abs.(instabilitiesUb[res,:] .- targetInstability))
        maxLambdas[res] = (lambdaRange)[xx[end]]  # to the right
        # note for typical bStARS, where you know what instability cutoff you
        # want you'd use the upperbound to find the min lambda and the lb to
        # find the max lambda    
    end

    netInstabilitiesUb = netInstabilitiesUb ./ totEdges
    netInstabilitiesLb = netInstabilitiesLb ./ totEdges
    maxLb = findmax(netInstabilitiesLb)
    maxLbInd = findall(x -> x == maxLb[1], netInstabilitiesLb)
    netInstabilitiesLb[maxLbInd[end]:end] .= maxLb[1] # take supremum for lambdas smaller than instability max
    maxUb = findmax(netInstabilitiesUb)
    maxUbInd = (findall(x -> x == maxUb[1], netInstabilitiesUb))
    netInstabilitiesUb[maxUbInd[end]:end] .= maxUb[1]
    xx = findmin(abs.(netInstabilitiesLb .- targetInstability))
    maxInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesLb .- targetInstability))
    minLambdaNet = (lambdaRange)[maxInstInd[end]]
    xx = findmin(abs.(netInstabilitiesUb .- targetInstability))
    minInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesUb .- targetInstability))
    maxLambdaNet = (lambdaRange)[minInstInd[end]]

    grnData.minLambdaNet = minLambdaNet
    grnData.maxLambdaNet = maxLambdaNet
    grnData.maxLambdas = maxLambdas
    grnData.minLambdas = minLambdas
    grnData.netInstabilitiesLb = netInstabilitiesLb
    grnData.netInstabilitiesUb = netInstabilitiesUb
    grnData.instabilitiesLb = instabilitiesLb
    grnData.instabilitiesUb = instabilitiesUb
    grnData.lambdaRangeWarm = lambdaRange   # ADDED: store warm-start λ grid for instability bound plots
end

function bstartsEstimateInstability(grnData::GrnData; totLambdas = 10, instabilityLevel = "Gene", zTarget = false, targetInstability::Float64 = 0.05, outputDir::Union{String, Nothing}=nothing)  # ADDED: targetInstability for λ selection diagnostic plot

    totResponses,totSamps = size(grnData.responseMat) # totResponses is same as length(grnData["targGenes"])
    totPreds = size(grnData.predictorMat,1)

    # if instabilityLevel == "Gene"
    #     minLambda = minimum(grnData.minLambdas)
    #     maxLambda = maximum(grnData.maxLambdas)
    # elseif instabilityLevel == "Network"
    #     minLambda = grnData.minLambdaNet
    #     maxLambda = grnData.maxLambdaNet
    # else
    #     println("Use either 'Gene' or 'Network' for instabilityLevel")
    # end

    # Gene Level Instabilities
    lambdaRangeGene = Vector{Vector{Float64}}(undef, totResponses)
    for res in 1:totResponses
        λmin = grnData.minLambdas[res]
        λmax = grnData.maxLambdas[res]
        lambdaRangeGene[res] = reverse(collect(range(λmin, stop=λmax, length=totLambdas)))
    end
    grnData.lambdaRangeGene = lambdaRangeGene

    # Network Level Instability
    minLambda = grnData.minLambdaNet
    maxLambda = grnData.maxLambdaNet

    lambdaRange = collect(range(minLambda, stop = maxLambda, length = totLambdas))
    lambdaRange = reverse(lambdaRange)
    grnData.lambdaRange = lambdaRange

    geneInstabilities = zeros(totResponses,totLambdas)
    netInstabilities = zeros(totLambdas,1)
    totEdges = 0   # denominator for network Instabilities
    # store number of subsamples for which an edge was nonzero, given that some
    # prior weights can be set to infinity, track to make sure these edges are
    # not counted
    ssMatrix = Inf*ones(totLambdas,totResponses,totPreds)
    subsamps = grnData.subsamps
    totSS = size(subsamps)[1]

    # get (finite) predictor indices for each response
    responsePredInds = Vector{Vector{Int}}(undef,0)
    for res = 1:totResponses
        currWeights = grnData.penaltyMat[res,:]
        push!(responsePredInds,findall(x -> x!=Inf, currWeights))
    end    

    totEdges = zeros(totResponses)
    betas = Array{Float64,3}(undef, totResponses, totPreds, totLambdas)
    Threads.@threads for res in ProgressBar(1:totResponses)
        lambdaRange = instabilityLevel == "Network" ?
                    grnData.lambdaRange :           # single shared vector
                    grnData.lambdaRangeGene[res]   # per‑gene vector
        predInds = responsePredInds[res]
        currPredNum = length(predInds)
        totEdges[res] = currPredNum 
        penaltyFactor = grnData.penaltyMat[res,predInds] 
        ssVals = zeros(totLambdas,currPredNum)
        for ss = 1:totSS
            subsamp = subsamps[ss,:]
            dt = fit(ZScoreTransform, grnData.predictorMat[predInds, subsamp], dims=2)
            currPreds = transpose(StatsBase.transform(dt, grnData.predictorMat[predInds, subsamp]))
            if zTarget
                dt = fit(ZScoreTransform, grnData.responseMat[res, subsamp], dims=1)
                currResponses = StatsBase.transform(dt, grnData.responseMat[res, subsamp])
            else
                currResponses = grnData.responseMat[res, subsamp]
            end
            lsoln = glmnet(currPreds, currResponses, penalty_factor = penaltyFactor, lambda = lambdaRange, alpha = 1.0)        
            currBetas = lsoln.betas
            betas[res,predInds, :] = currBetas
            ssVals = ssVals + abs.(sign.(currBetas))' 
        end
        ssMatrix[:,res,predInds] = ssVals
    end

    for res = 1:totResponses    
        currWeights = grnData.penaltyMat[res,:]
        predInds = responsePredInds[res]
        currPredNum = length(predInds)
        ssVals = zeros(totLambdas,currPredNum)
        ssVals[:,:] = ssMatrix[:,res,predInds]   
        theta2 = (1/totSS)*ssVals # empirical edge probability, lambdas X currPreds
        instabilitiesPerEdge = 2*(theta2 .* (1 .- theta2))
        aveInstabilities = mean(instabilitiesPerEdge,dims=2) # lambdas X 1    
        maxUb = (findmax(aveInstabilities))
        instabSUP = aveInstabilities
        maxUbInd = first(Tuple(maxUb[2]))
        instabSUP[maxUbInd:end,1] .= maxUb[1]
        geneInstabilities[res,:] = instabSUP
    end
    
    ## calculate instabilities network-wise
    currSS = zeros(totResponses,totPreds)
    instabRange = zeros(totLambdas,1)
    for lind = 1:totLambdas # start at highest lambda (lowest instability and work down)
        currSS[:,:] = ssMatrix[lind,:,:]    
        theta2 = (1/totSS)*currSS # empirical edge probability, responses X currPreds
        instabilitiesPerEdge = 2*(theta2 .* (1 .- theta2))
        instabVec = instabilitiesPerEdge[:]
        validEdges = findall(isfinite.(currSS[:])) # limit to finite edges
        instabMax = findmax(instabRange)[1]
        netInstabilities[lind] = max(mean(instabVec[validEdges]),max(instabMax))
    end
    grnData.netInstabilities = vec(netInstabilities)
    grnData.geneInstabilities = geneInstabilities
    grnData.stabilityMat = ssMatrix
    grnData.betas = betas

    if outputDir !== nothing && outputDir !== ""
        outputFile = joinpath(outputDir, "instabOutMat.jld")
        save_object(outputFile, grnData)
        # ADDED: auto-generate network-level diagnostic plots alongside saved data
        plotInstabilityCurves(grnData;
                              mode              = :network,
                              targetInstability = targetInstability,
                              outputDir         = outputDir)
    end
end
