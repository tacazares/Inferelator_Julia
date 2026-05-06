# mutable struct BuildGrn
#         networkStability::Matrix{Float64}
#         lambda::Union{Float64, Vector{Float64}}
#         targs::Vector{String}
#         regs::Vector{String}
#         rankings::Vector{Float64}
#         signedQuantile::Vector{Float64}
#         partialCorrelation::Vector{Float64}
#         inPrior::Vector{String}
#         networkMat::Matrix{Any}
#         networkMatSubset::Matrix{Any}
#         inPriorVec::Vector{Float64}
#         betas::Matrix{Float64}
#         function BuildGrn()
#             return new(
#                 Matrix{Float64}(undef, 0, 0),   # networkStability
#                 0.0,                            # lambda
#                 [],                             # targs
#                 [],                             # regs
#                 [],                             # rankings
#                 [],                             # signedQuantile
#                 [],                             # partialCorrelation 
#                 [],                             # inPrior
#                 Matrix{Float64}(undef, 0, 0),   # networkMat
#                 Matrix{Float64}(undef, 0, 0),   # networkMatSubset
#                 [],                             # inPriorVec
#                 Matrix{Float64}(undef, 0, 0)    # betas
#             )
#         end
#     end

function chooseLambda!(grnData::GrnData, buildGrn::BuildGrn; instabilityLevel = "Gene", targetInstability = 0.05)
    totLambdas, totNetGenes, totNetTfs = size(grnData.stabilityMat)
    networkStability = zeros(totNetGenes, totNetTfs)
    networkStability = convert(Matrix{Float64}, networkStability)
    # lambdaRange = reverse(grnData.lambdaRange)
    lambdaRange = grnData.lambdaRange
    betas = zeros(totNetGenes, totNetTfs)
    lambdaTrack = []
    ## transform StARS instabilities into stabilities 
    if instabilityLevel == "Gene"
        totMins = 0
        totMaxs = 0
        for targ = 1:totNetGenes
            currInstabs = grnData.geneInstabilities[targ, :]
            devs = abs.(currInstabs .- targetInstability)
            globalMin = minimum(devs)
            minInds = findall(x -> x == globalMin[1], devs) # globalMin
            minInd = minInds[end]

            lambdaRangeGene = grnData.lambdaRangeGene[targ]
            push!(lambdaTrack, lambdaRangeGene[minInd])
        
            networkStability[targ, :] = grnData.stabilityMat[minInd, targ, :]
            betas[targ, :] = grnData.betas[targ, :, minInd]

            if minInd == 1
                totMins += 1
            elseif minInd == totLambdas
                totMaxs += 1
            end
        end
        if totMins > 0
            @info "Target instability reached at minimum lambda for $totMins gene(s)"
        end
        if totMaxs > 0
            @info "Target instability reached at maximum lambda for $totMaxs gene(s)"
        end

    elseif instabilityLevel == "Network"
        @info "Network instabilities detected"
        # find the single lambda corresponding to the cutoff
        devs = abs.(grnData.netInstabilities .- targetInstability)
        globalMin = findmin(devs)
        minInds = first(Tuple(globalMin[2]))
        globalMin = globalMin[1]
        # take the largest lambda that is closest to targetInstability
        minInd = minInds[end]
        if minInd == 1
            @info "Minimum lambda used" instability=grnData.netInstabilities[minInd] targetInstability
        elseif minInd == totLambdas
            @info "Maximum lambda used" instability=grnData.netInstabilities[minInd] targetInstability
        end
        networkStability[:, :] = grnData.stabilityMat[minInd, :, :]
        betas = grnData.betas[:, :, minInd]
        lambdaTrack = lambdaRange[minInds]

    else
        error("instabSource not recognized, should be either Gene or Network.")
    end

    # ssMatrix has infinity entries to mark illegal TF-gene interactions
    # (e.g., TF mRNA TFA cannot be used to predict TF gene expression)
    networkStabilityVec = networkStability[:]
    networkStabilityVec[isinf.(networkStabilityVec)] .= 0
    networkStability = reshape(networkStabilityVec, totNetGenes, totNetTfs)

    if length(lambdaTrack) == 1
        lambdaTrack = convert(Float64, lambdaTrack)
    else
        lambdaTrack = convert(Vector{Float64}, lambdaTrack)
    end

    buildGrn.lambda = lambdaTrack
    buildGrn.networkStability = networkStability
    buildGrn.betas = betas
    @info "chooseLambda! selected" chosenλ = lambdaTrack instabilityLevel = instabilityLevel
end

# function firstNByGroup(vect::AbstractVector, N::Integer)
#     """
#     firstNByGroupIndices(vec::AbstractVector, N::Integer)

#     Return the indices of the first `N` occurrences of each unique element in `vec`.

#     This is useful when you want to subset another array based on limited occurrences per group.

#     # Arguments
#     - `vec`: A vector of values to group by (e.g., transcription factors).
#     - `N`: The maximum number of elements to select per unique group.

#     # Returns
#     - A vector of indices corresponding to the first `N` entries per group.
#     """
#     seen = Dict{eltype(vect), Int}()
#     idxs = Int[]
#     for (i, v) in enumerate(vect)
#         seen[v] = get(seen, v, 0) + 1
#         if seen[v] <= N
#             push!(idxs, i)
#         end
#     end
#     return idxs
# end

function rankEdges!(expressionData::GeneExpressionData, tfaData::PriorTFAData, grnData::GrnData, buildGrn::BuildGrn; 
                    mergedTFsData::Union{mergedTFsResult, Nothing}=nothing,
                    useMeanEdgesPerGeneMode = true, 
                    meanEdgesPerGene = 20, 
                    correlationWeight = 1,
                    outputDir::Union{String, Nothing}=nothing)
    
    if isempty(outputDir)
        error("You must provide a non-empty outputDir to save the results.")
    end

    predictorMat = grnData.predictorMat
    allPredictors = grnData.allPredictors

    # stabilityMat is populated by bStARS but not by EBIC/bEBIC — fall back to
    # networkStability dimensions when it is empty.
    if !isempty(grnData.stabilityMat) && ndims(grnData.stabilityMat) == 3
        totLambdas, totNetGenes, totNetTfs = size(grnData.stabilityMat)
    else
        totNetGenes, totNetTfs = size(buildGrn.networkStability)
        totLambdas = 0
    end
    totInts = totNetGenes * totNetTfs 
    targs = repeat(expressionData.targGenes, totNetTfs, 1)
    regs = vec(repeat(permutedims(allPredictors), totNetGenes,1))
    rankTmp = buildGrn.networkStability[:]
    keepInds = findall(x -> x != 0 && x != Inf, rankTmp) # keep nonzero, remove infinite values (e.g., corresponding to TF-TF edges when TF mRNA used for TFA)
    rankTmp = rankTmp[keepInds]
    inds = reverse(sortperm(rankTmp))
    rankings = sort(rankTmp,rev=true)
    regs = regs[keepInds[inds]]
    targs = targs[keepInds[inds]]
    totInfInts = length(rankings)

    if useMeanEdgesPerGeneMode
        totQuantEdges = length(unique(targs)) * meanEdgesPerGene
        selectionIndices = 1:min(totQuantEdges, totInfInts)
    else
        selectionIndices = firstNByGroup(targs, meanEdgesPerGene)
    end

    # Compute quantiles 
    F = ecdf(rankings[selectionIndices])  # uses the StatsBase package
    quantiles = F.(rankings[selectionIndices])

    ## take what's in the meanEdgesPerGene network and get partial correlations
    allCoefs = zeros(totNetGenes,totNetTfs)
    allQuants = zeros(totNetGenes,totNetTfs)
    allStabsTest = buildGrn.networkStability
    keptTargs = permutedims(targs[selectionIndices])
    uniTargs = unique(keptTargs)
    totUniTargs = length(uniTargs)
    tfsPerGene = zeros(totUniTargs,1)

    for targ = ProgressBar(1:totUniTargs)
        currTarg = uniTargs[targ]
        targRankInds = last.(Tuple.(findall(x -> x == currTarg, keptTargs)))
        currRegs = regs[targRankInds]

        # --- Deduplicate currRegs while keeping aligned targRankInds ---
        uniCurrRegsInds = firstNByGroup(currRegs, 1)  # This reurns the index where each unique currRegs first appeared
        uniCurrRegs = currRegs[uniCurrRegsInds]    
        targRankInds = targRankInds[uniCurrRegsInds]  # Ensure the target corresponding to the duplicated Regulator is removed
        # --- Get target index in full gene list
        targInd = last.(Tuple.(findall(x -> x == currTarg, expressionData.targGenes)))
        tfsPerGene[targ] = Int.(length(targRankInds))
        # tfsPerGene = Int.(tfsPerGene)

        # --- Match predictors ---
        matchedRegs = intersect(allPredictors,currRegs)
        # inds = findall(in(matchedRegs), allPredictors)
        # regressIndsMat = first.(inds)
        regressIndsMat = findall(in(matchedRegs), allPredictors)
        rankVecInds = findall(in(matchedRegs),uniCurrRegs)
        # --- Prepare input matrix: target + predictors ---
        currTargVals = vec(transpose(grnData.responseMat[targInd,:]))
        currPredVals = transpose(predictorMat[regressIndsMat,:])
        combTargPreds = vcat(currTargVals', currPredVals')'
        combTargPreds = permutedims(combTargPreds')
        prho = partialCorrelationMat(combTargPreds; first_vs_all = true)
        prho = prho[2:end]
        prho = vec(prho)

        if length(findall(x -> x == NaN, prho)) == 0  # make sure there weren't too many edges, 
            allCoefs[targInd,regressIndsMat] = prho
        else
            @warn "Partial correlation was singular for $currTarg" nTFs=length(regressIndsMat)
        end   

        allQuants[targInd,regressIndsMat] = quantiles[targRankInds[rankVecInds]]
        allStabsTest[targInd,regressIndsMat] = allStabsTest[targInd,regressIndsMat] + (correlationWeight .* transpose(round.(abs.(prho),digits = 4)))
    end

    inPriorMat = sign.(abs.(grnData.priorMatProcessed))
    
    mergeTfLocVec = zeros(totNetTfs) # for keeping track of merged TFs (needed for partial correlation calculation)
    if (mergedTFsData !== nothing) && 
        (mergedTFsData.mergedPrior !== nothing) && 
        (mergedTFsData.mergedTFMap !== nothing)
        
        mergedTFs =  mergedTFsData.mergedTFMap[:, 1]
        individualTFs =  mergedTFsData.mergedTFMap[:, 2]
        totMerged = length(individualTFs)

        rmInds = []        # remove merged TFs from regulators
        addRegs = []
        addInts = []
        addCoefs = []
        addPMat = []
        addPredMat = []
        addLoc = []
        addQuants = []
        for mind = 1:totMerged
            mTf = mergedTFs[mind]
            inputLocs = findall(x -> x == mTf,allPredictors)
            totMInts = length(inputLocs) # number of interactions for merged TF
            rmInds = [rmInds; inputLocs];
            if length(inputLocs) > 0
                indTfs = intersect(permutedims(split(individualTFs[mind],", ")),tfaData.pRegsNoTfa); # intersect ensures that TF was a potential regulator (e.g., based on gene expression)
                totIndTfs = length(indTfs)
                for indt = 1:totIndTfs
                    indTf = indTfs[indt]
                    append!(addRegs, fill(indTf, 1))
                    append!(addInts, allStabsTest[:,inputLocs])
                    append!(addPMat, inPriorMat[:,inputLocs])
                    append!(addQuants, allQuants[:,inputLocs])
                    append!(addCoefs, allQuants[:,inputLocs])
                    append!(addPredMat, predictorMat[inputLocs,:])
                    append!(addLoc, mind)
                end                
            end
        end
        @info "$(length(rmInds)) merged TFs expanded to individual TFs"
        keepInds = setdiff(1:totNetTfs,rmInds)
        # remove merged TFs and add individual TFs
        if length(addRegs) > 0
            allPredictors = collect(allPredictors[keepInds])
            append!(allPredictors, addRegs)
            allStabsTest = hcat(allStabsTest[:,keepInds], reshape(addInts, size(allStabsTest)[1], Int(size(addInts)[1] / size(allStabsTest)[1])))
            allCoefs = hcat(allCoefs[:,keepInds], reshape(addCoefs, size(allStabsTest)[1], Int(size(addInts)[1] / size(allStabsTest)[1])))
            allQuants = hcat(allQuants[:,keepInds], reshape(addQuants, size(allStabsTest)[1], Int(size(addInts)[1] / size(allStabsTest)[1])))
            inPriorMat = hcat(inPriorMat[:,keepInds], reshape(addPMat, size(inPriorMat)[1], Int(size(addPMat)[1] / size(inPriorMat)[1])))
            predictorMat = vcat(predictorMat[keepInds, :], reshape(addPredMat, Int((size(addPredMat)[1]) / size(predictorMat)[2]), size(predictorMat)[2]))
            append!(mergeTfLocVec, addLoc)
        end
    else
        @info "No merged TFs data found — skipping TF expansion"
    end

    ## re-rank based on possibly de-merged TFs
    rankings = allStabsTest[:]       
    coefVec = allCoefs[:]

    
    inPriorVec = inPriorMat[:]
    totNetTfs = length(allPredictors)
    totInts = totNetGenes * totNetTfs
    targs = repeat((expressionData.targGenes),totNetTfs,1)
    regs1 = repeat(permutedims(allPredictors),totNetGenes,1)
    regs = reshape(regs1,totInts,1)

    rankings = convert(Vector{Float64}, rankings)
    coefVec = convert(Vector{Float64}, coefVec)
    inPriorVec = convert(Vector{Float64}, inPriorVec)
    predictorMat = convert(Matrix{Float64}, predictorMat)
    allStabsTest = convert(Matrix{Float64}, allStabsTest)
    allCoefs = convert(Matrix{Float64}, allCoefs)
    allQuants = convert(Matrix{Float64}, allQuants)
    inPriorMat = convert(Matrix{Float64}, inPriorMat)

    # Sync back — allStabsTest may have been rebound by hcat (merged TF expansion)
    # and convert, so buildGrn.networkStability must be updated explicitly here.
    buildGrn.networkStability = allStabsTest

    ## only keep nonzero rankings
    keepRankings = findall(x -> x != 0 && x != Inf, rankings)
    indsMerged = sortperm(rankings[keepRankings])
    indsMerged = reverse(indsMerged)

    # update info sources
    rankings = rankings[keepRankings[indsMerged]]
    coefVec = coefVec[keepRankings[indsMerged]]
    inPriorVec = inPriorVec[keepRankings[indsMerged]]    
    regs = regs[keepRankings[indsMerged]]
    targs = targs[keepRankings[indsMerged]]
    totInfInts = length(rankings)

    # totQuantEdges = length(unique(targs))*meanEdgesPerGene
    # Compute quantiles
    F = ecdf(rankings)  
    quantiles = F.(rankings)

    # Compute signedQuantiles
    signedQuantile = sign.(coefVec) .* quantiles

    ## ------- Color calculations for JP-Gene-Viz
    # Compute strokeWidth and colors
    if isempty(rankings)
        error("rankEdges!: no non-zero edges remain after filtering — the inferred " *
              "network is empty. Likely cause: model selection chose an overly sparse " *
              "solution (all coefficients zero). " *
              "Common triggers for EBIC / bEBIC: " *
              "(1) zScoreLASSO=false without explicit minLambda/maxLambda — the response " *
              "is on its raw scale so lambda_max is inflated and EBIC picks a null model; " *
              "(2) gamma too high (current default 0.5) — try gamma=0.0 or gamma=0.25; " *
              "(3) lambda range too narrow — pass minLambda/maxLambda or leave as nothing " *
              "to let GLMNet choose. " *
              "See ebicSelect! / bebicEstimateLambdas! docstrings for details.")
    end
    minRank, maxRank = extrema(rankings)
    rankRange = max(maxRank - minRank, eps())  # prevent division by zero
    # Dash Styling
    strokeDashArray = ifelse.(inPriorVec .!= 0, "Yes", "No")
    networkMatrix = hcat(regs, targs, signedQuantile, rankings, coefVec, strokeDashArray)
    if useMeanEdgesPerGeneMode
        totQuantEdges = length(unique(targs)) * meanEdgesPerGene
        # selectionIndices = 1:totQuantEdges
        selectionIndices = 1:min(totQuantEdges, totInfInts)
    else
        selectionIndices = firstNByGroup(targs, meanEdgesPerGene)
        
    end
    
    networkMatrixSubset = hcat(regs[selectionIndices], targs[selectionIndices], signedQuantile[selectionIndices],
                                rankings[selectionIndices], coefVec[selectionIndices], strokeDashArray[selectionIndices]
                            )
    
    buildGrn.regs = regs
    buildGrn.targs = targs
    buildGrn.signedQuantile = signedQuantile
    buildGrn.rankings = rankings
    buildGrn.partialCorrelation = coefVec
    buildGrn.inPrior = strokeDashArray
    buildGrn.networkMat = networkMatrix
    buildGrn.networkMatSubset = networkMatrixSubset
    buildGrn.inPriorVec = inPriorVec
    buildGrn.mergeTfLocVec = mergeTfLocVec   # just added


    if outputDir !== nothing && outputDir !== ""
        outputFile = joinpath(outputDir, "grnOutMat.jld")
        save_object(outputFile, buildGrn)
    end
end

# function saveData(expressionData::GeneExpressionData, tfaData::PriorTFAData, grnData::GrnData, buildGrn::BuildGrn, outputDir::String, fileName::String)
#     output = outputDir * "/" * fileName
#     @save output expressionData tfaData grnData buildGrn
# end

# function writeNetworkTable!(buildGrn::BuildGrn; outputDir::String, networkName::Union{String, Nothing}=nothing)
#     local baseName = (networkName === nothing || isempty(networkName)) ? "edges" : networkName * "edges"

#     outputFile = joinpath(outputDir, baseName * ".tsv")
#     colNames = "TF\tGene\tsignedQuantile\tStability\tCorrelation\tinPrior\n"
#     open(outputFile, "w") do io
#         write(io, colNames)
#         writedlm(io, buildGrn.networkMat)
#     end

#     outputFileSubset = joinpath(outputDir, baseName * "_subset.tsv")
#     open(outputFileSubset, "w") do io
#         write(io, colNames)
#         writedlm(io, buildGrn.networkMatSubset)
#     end
# end
