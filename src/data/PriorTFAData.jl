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

            @info "Case 1: TF mRNA — using unmerged prior, limiting to potential regulators" 
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
            # output mRNA estimates for regulators without TFA
            noPriorRegs = setdiff(expressionData.potRegsmRNA, pRegsNoTfa)
            expInds = findall(in(noPriorRegs), expressionData.potRegsmRNA)
            noPriorRegsMat = expressionData.potRegMatmRNA[expInds,:]

            @info "Case 2: prior-based TFA" 
            ## Case 2: prior-based TFA is used
            # Check whether there were degenerate TFs and outputs
            if (mergedTFsData !== nothing) &&
                (mergedTFsData.mergedPrior !== nothing) &&
                (mergedTFsData.mergedTFMap !== nothing)

                @info "Using merged degenerate TFs prior file" 
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
            # remove predictors that have fewer than the specified minimum number of targets
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
                @warn "Gene order in expression matrix and prior matrix do not match" 
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
            @info "Target expression normalized (z-score per gene)" 
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
                    # BUG (fixed): rand(targInds, ssampleSize) samples WITH replacement,
                    # so the same index can appear twice, silently overwriting rows in sPrior
                    # and reducing the effective number of unique edges below the intended ~63%.
                    # MATLAB's randperm(n, k) always samples WITHOUT replacement.
                    # Fix: StatsBase.sample(..., replace=false)
                    ssample = StatsBase.sample(targInds, ssampleSize, replace=false)

                    sPrior[ssample, col] = priorMatrix[ssample, col]
                end
                tfas[ss, :, :] = sPrior \ targExp
            end

            # BUG (fixed): median(tfas, dims=1) returns a 3D array of shape (1, totPreds, totConds)
            # because Julia keeps the reduced dimension as a singleton. Assigning that to
            # medTfas (typed Matrix{Float64}) either errors or stores the wrong shape.
            # MATLAB's median(tfas, 1) followed by (:,:) indexing squeezes the singleton dim.
            # Fix: dropdims(..., dims=1) to get the expected (totPreds, totConds) matrix.
            priorData.medTfas = dropdims(median(tfas, dims=1), dims=1)
            @info "Median from $edgeSS subsamples used for prior-based TFA" 
        else
            # solves for X = Prior * TFA.
            # TFA = argmin||priorMatrix * TFA - targExp||²
            priorData.medTfas = priorMatrix \ targExp
            @info "No subsampling — all edges used for prior-based TFA estimate" 
        end

        if outputDir !== nothing && outputDir !== ""

            outputFile = joinpath(outputDir, "tfaMat.jld")
            save_object(outputFile, priorData)

        end
    end


    """
        applyTimeLag!(priorData, expressionData, timeLagFile, timeLag)

    Apply a time-lag correction to TF activity and mRNA estimates after TFA has
    been calculated.

    For each consecutive time-point pair (P → Q) listed in `timeLagFile`, the
    correction adjusts the value at sample Q by subtracting the extrapolated
    change due to the lag:

        adjusted[Q] = original[Q] - (original[Q] - original[P]) / Δt * timeLag

    where Δt = timeQ − timeP (in the same units as `timeLag`).

    Three matrices are corrected in-place:
    - `priorData.medTfas`           — prior-based TFA estimates
    - `priorData.noPriorRegsMat`    — mRNA for regulators absent from the prior
    - `expressionData.potRegMatmRNA` — mRNA for all potential regulators

    # Arguments
    - `priorData::PriorTFAData`         : TFA struct (must have `medTfas` populated)
    - `expressionData::GeneExpressionData` : expression struct with `cellLabels`
    - `timeLagFile::String`             : path to a 4-column TSV with a header row.
      Columns: SampleQ name, time of Q, SampleP name, time of P.
      Only samples that appear as time-series pairs need to be listed — the rest
      are left unchanged.
    - `timeLag::Real`                   : expected lag between TF mRNA and protein
      activity, in the same units as the time values in `timeLagFile`.

    # Notes
    - Based on Bonneau et al. (2006) Genome Biology and the MATLAB function
      `integratePrior_estTFA_timeLag` from the Miraldi lab.
    - Rows in `timeLagFile` that name samples not found in `cellLabels` are
      skipped with a warning rather than erroring out.
    """
    function applyTimeLag!(priorData::PriorTFAData,
                           expressionData::GeneExpressionData,
                           timeLagFile::String,
                           timeLag::Real)

        isempty(timeLagFile) && return

        # ── 1. Parse the time-lag file ────────────────────────────────────────
        # Expected format (tab-separated, one header line):
        #   SampleQ  TimeQ  SampleP  TimeP
        #   where SampleP immediately precedes SampleQ in the time series.
        lines = readlines(timeLagFile)
        timePointConds        = String[]
        timePointProceedConds = String[]
        distanceBetweenPoints = Float64[]

        for line in lines[2:end]           # skip header
            parts = split(line, '\t')
            length(parts) < 4 && continue
            push!(timePointConds,        strip(parts[1]))
            timeQ = parse(Float64, strip(parts[2]))
            push!(timePointProceedConds, strip(parts[3]))
            timeP = parse(Float64, strip(parts[4]))
            push!(distanceBetweenPoints, timeQ - timeP)
        end
        totTimes = length(timePointConds)

        # ── 2. Snapshot unmodified matrices (slopes use the pre-lag values) ──
        # The MATLAB uses *Tmp copies for slope calculation while writing
        # corrections into the originals, so all pairs see the same baseline.
        noPriorRegsMatOrig   = copy(priorData.noPriorRegsMat)
        medTfasOrig          = copy(priorData.medTfas)
        potRegMatmRNAOrig    = copy(expressionData.potRegMatmRNA)

        # ── 3. Apply correction for each time-point pair ─────────────────────
        for tp in 1:totTimes
            currInd    = findfirst(==(timePointConds[tp]),        expressionData.cellLabels)
            proceedInd = findfirst(==(timePointProceedConds[tp]), expressionData.cellLabels)

            if isnothing(currInd) || isnothing(proceedInd)
                @warn "applyTimeLag!: sample(s) not found in cellLabels — skipping pair " *
                      "($(timePointConds[tp]), $(timePointProceedConds[tp]))"
                continue
            end

            dist = distanceBetweenPoints[tp]

            # TF mRNA for regulators not in the prior
            delta = noPriorRegsMatOrig[:, currInd] .- noPriorRegsMatOrig[:, proceedInd]
            priorData.noPriorRegsMat[:, currInd] .-= (delta ./ dist) .* timeLag

            # Prior-based TFA
            delta = medTfasOrig[:, currInd] .- medTfasOrig[:, proceedInd]
            priorData.medTfas[:, currInd] .-= (delta ./ dist) .* timeLag

            # All TF mRNA
            delta = potRegMatmRNAOrig[:, currInd] .- potRegMatmRNAOrig[:, proceedInd]
            expressionData.potRegMatmRNA[:, currInd] .-= (delta ./ dist) .* timeLag
        end

        @info "Time-lag correction applied to $totTimes time-point pairs (lag = $timeLag)" 
    end
