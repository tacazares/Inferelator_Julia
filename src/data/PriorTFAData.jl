    # Standard libraries
    using LinearAlgebra
    using Statistics
    using DelimitedFiles
    using JLD2

    # Struct defined in src/Types.jl

    """
        processPriorFile!(priorData::PriorTFAData, priorFile)

    Reads and processes a prior file, storing the extracted information in a `PriorTFAData` object.

    # Arguments
    - `priorData::PriorTFAData`: The struct to be populated with data from the prior file.
    - `priorFile::String`: The path to the prior file, formatted with tab-separated values.

    # Updates
    - `priorData.pRegs`: Stores the list of transcription factors (TFs) from the prior file.
    - `priorData.pTargs`: Stores the list of target genes from the prior file.
    - `priorData.priorMatrix`: Stores the interactions matrix, indicating TF-gene relationships.

    # Raises
    - An error if the prior file does not exist.

    # Notes
    - Assumes the first row in the file contains TF names, and subsequent rows represent interactions with target genes.
    - The gene list and matrix are sorted alphabetically by gene names.
    """
    function processPriorFile!(priorData::PriorTFAData,
                            expressionData::GeneExpressionData,
                            priorFile;  mergedTFsData::Union{mergedTFsResult, Nothing}=nothing, minTargets = 3)

        if isfile(priorFile)

            println("--- Case 1")
            fid = open(priorFile)
            C = readdlm(fid, '\t', '\n', skipstart=0)
            close(fid)

            # Process and store genes and interactions matrix
            pRegsTmp = convert(Vector{String}, filter(!isempty, C[1, :]))
            C = C[2:end, :]
            inds = sortperm(C[:, 1])
            C = C[inds, :]
            pTargsTmp = convert(Vector{String}, C[:, 1])
            pMatrixTmp = convert(Matrix{Float64}, C[:, 2:end])

            # Filter to only include those in potential regulator and target gene list
            targInds = findall(in(expressionData.tfaGenes), pTargsTmp)
            regInds = findall(in(expressionData.potRegs), pRegsTmp)
            pRegsNoTfa = pRegsTmp[regInds]
            pTargsNoTfa = pTargsTmp[targInds]
            priorMatrixNoTfa = pMatrixTmp[targInds, regInds]

            # Find TFs that have expression data but arent in the prior
            noPriorRegs = setdiff(expressionData.potRegsmRNA, pRegsNoTfa)
            expInds = findall(in(noPriorRegs), expressionData.potRegsmRNA)
            noPriorRegsMat = expressionData.potRegMatmRNA[expInds,:]

            println("--- Case 2")
            ## Case 2: prior-based TFA is used
            # Check whether there were degenerate TFs and outputs
            if (mergedTFsData !== nothing) &&
                (mergedTFsData.mergedPrior !== nothing) &&
                (mergedTFsData.mergedTFMap !== nothing)

                println("-------- Using merge degenerate TFs prior file")
                mergedTFs =  mergedTFsData.mergedTFMap[:, 1]
                individualTFs =  mergedTFsData.mergedTFMap[:, 2]

                totMergedSets = length(individualTFs)
                keepMergedIndices = Int[]

                for idx in 1:totMergedSets
                    currSet = split(individualTFs[idx], ", ")
                    usedTfs = intersect(currSet, expressionData.potRegs)
                    if !isempty(usedTfs)
                        push!(keepMergedIndices, idx)
                    end
                end

                if !isempty(keepMergedIndices)  # add merged potential regulators to our list
                    expressionData.potRegs = union(expressionData.potRegs, mergedTFs[keepMergedIndices])
                    # Now load the merged prior matrix data (mergedPrior)
                    priorDF = mergedTFsData.mergedPrior
                    rowInd = sortperm(priorDF[:, 1])
                    pRegsTmp = names(priorDF)[2:end]
                    totPRegs = length(pRegsTmp)
                    pTargsTmp = priorDF[rowInd, 1]
                    pMatrixTmp = Matrix(priorDF[rowInd, 2:end])
                end
            end

            ### Filter TFs and Genes for TFA calculation
            pTargInds = findall(in(expressionData.tfaGenes), pTargsTmp)
            pRegInds = findall(in(expressionData.potRegs), pRegsTmp)
            pTargs = pTargsTmp[pTargInds]
            pRegs = pRegsTmp[pRegInds]
            pInts = pMatrixTmp[pTargInds,pRegInds]
            # Filter for minimum target genes in prior
            interactionsPerTF = sum((abs.(sign.(pInts))), dims=1)
            keepRegs = Tuple.(findall(x -> x > minTargets, interactionsPerTF))
            keepRegs = last.(keepRegs)
            pInts = pInts[:,keepRegs]
            pRegs = pRegs[keepRegs]
            # Filter genes with no regulators
            interactionsPerTarg = sum(abs.(sign.(pInts)), dims = 2)
            keepTargs = Tuple.(findall(x -> x > 0, interactionsPerTarg))
            keepTargs = first.(keepTargs)
            pInts = pInts[keepTargs,:]
            pTargs = pTargs[keepTargs]

            ### Ensure expression and prior targets are in the same order
            expTargInds = findall(in(pTargs), expressionData.tfaGenes)
            targExp = expressionData.tfaGeneMat[expTargInds,:]
            if !(pTargs == expressionData.tfaGenes[expTargInds])
                println("Warnings, gene order in expression matrix and prior matrix do not match!!")
            end

            # Add data to object
            priorData.pRegs = pRegs
            priorData.pTargs = pTargs
            priorData.priorMatrix = pInts
            priorData.pRegsNoTfa = pRegsNoTfa
            priorData.pTargsNoTfa = pTargsNoTfa
            priorData.priorMatrixNoTfa = priorMatrixNoTfa
            priorData.noPriorRegs = noPriorRegs
            priorData.noPriorRegsMat = noPriorRegsMat
            priorData.targExpression = targExp

        else
            error("Prior file not found.")
        end
    end

    function calculateTFA!(priorData::PriorTFAData, expressionData::GeneExpressionData;
                            edgeSS = 0, zTarget::Bool = false, outputDir::Union{String, Nothing}=nothing)
        priorMatrix = priorData.priorMatrix
        targExp = priorData.targExpression
        totTargs = size(priorMatrix, 1)
        totPreds = size(priorMatrix, 2)
        totConds = size(targExp, 2)

        if zTarget
            targExp = (targExp .- mean(targExp, dims=2)) ./ std(targExp, dims=2)
            println("Target expression normalized (z-score per gene).")
        end

        if edgeSS > 0
            tfas = zeros(Float64, edgeSS, totPreds, totConds)

            for ss in 1:edgeSS
                sPrior = zeros(Float64, totTargs, totPreds)

                for col in 1:totPreds
                    currTargs = priorMatrix[:, col]
                    targInds = findall(!iszero, currTargs)
                    totCurrTargs = length(targInds)
                    ssampleSize = Int(ceil(0.63 * totCurrTargs))
                    ssample = rand(targInds, ssampleSize)

                    sPrior[ssample, col] = priorMatrix[ssample, col]
                end
                tfas[ss, :, :] = sPrior \ targExp
            end

            priorData.medTfas = median(tfas, dims=1)
            println("Median from ", string(edgeSS), " subsamples used for prior-based TFA.")
        else
            # solves for X = Prior * TFA.
            # TFA = argmin||priorMatrix * TFA - targExp||²
            priorData.medTfas = priorMatrix \ targExp
            println("No subsampling for prior-based TFA estimate.")
        end

        if outputDir !== nothing && outputDir !== ""

            outputFile = joinpath(outputDir, "tfaMat.jld")
            save_object(outputFile, priorData)

        end
    end
