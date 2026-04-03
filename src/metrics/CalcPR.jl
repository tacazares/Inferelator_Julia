"""
    computeMacroMetrics(gsPrecisionsByTf, gsRecallsByTf, gsFprsByTf, gsAuprsByTf, gsArocsByTf;
                        target_points=1000, min_step=1e-4, step_method=:min_gap)

Compute macro-averaged PR and ROC curves across all TFs using interpolation.

Per-TF curves are interpolated onto a common grid and averaged. Only TFs with
non-empty, non-zero recall arrays contribute to the average.

# Arguments
- `gsPrecisionsByTf` : Vector of per-TF precision arrays.
- `gsRecallsByTf`    : Vector of per-TF recall arrays.
- `gsFprsByTf`       : Vector of per-TF false positive rate arrays.
- `gsAuprsByTf`      : Vector of per-TF AUPR values.
- `gsArocsByTf`      : Vector of per-TF AUROC values.

# Keyword Arguments
- `target_points::Int=1000`       : Number of interpolation points if `step_method=:target_points`.
- `min_step::Float64=1e-4`        : Minimum allowable step size for interpolation.
- `step_method::Symbol=:min_gap`  : Step selection method — `:min_gap` or `:target_points`.

# Returns
An `OrderedDict` with two nested `OrderedDict`s:
- `:macroPR`  -> `:auprInterpolated`, `:precisions`, `:recalls`
- `:macroROC` -> `:aurocInterpolated`, `:fprs`, `:tprs`

# Example
```julia
results = computeMacroMetrics(precsByTf, recsByTf, fprsByTf, auprsByTf, aurocsByTf)
results[:macroPR][:auprInterpolated]   # macro AUPR
results[:macroPR][:precisions]         # macro precision curve
results[:macroROC][:tprs]              # macro TPR curve
```
"""
function computeMacroMetrics(gsPrecisionsByTf, gsRecallsByTf, gsFprsByTf, gsAuprsByTf, gsArocsByTf;
                               target_points=1000, 
                               min_step=1e-4,
                               step_method=:min_gap)

    # ---  Scalar Macro Metrics ---
    macroAUPR_scalar  = mean(filter(isfinite, gsAuprsByTf))
    macroAUROC_scalar = mean(filter(isfinite, gsArocsByTf))


    # --- Interpolated Macro Curves ---

    # Determine dynamic max recall/FPR
    allRecalls = reduce(vcat, filter(!isempty, gsRecallsByTf))
    allRecalls = filter(x -> isfinite(x), allRecalls) 
    allFprs = reduce(vcat, filter(!isempty, gsFprsByTf))
    allFprs = filter(x -> isfinite(x), allFprs)  
    
    maxRecall = isempty(allRecalls) ? 0.0 : maximum(allRecalls)
    maxFpr    = isempty(allFprs) ? 0.0 : maximum(allFprs)
    
    recStep = adaptiveStep(gsRecallsByTf; target_points=target_points, min_step=min_step, method=step_method)
    fprStep = adaptiveStep(gsFprsByTf;   target_points=target_points, min_step=min_step, method=step_method)
    recInterpPts = 0.0:recStep:maxRecall   # Interpolation points
    fprInterpPts = 0.0:fprStep:maxFpr

    macroPrecisions = zeros(length(recInterpPts))
    macroTprs       = zeros(length(fprInterpPts))
    nTFs = length(gsPrecisionsByTf)
    validTFs = [i for i in 1:nTFs if !isempty(gsRecallsByTf[i]) && maximum(gsRecallsByTf[i]) > 0]
    for i in validTFs
        # PR Curve
        interpPrec = dedupNInterpolate(gsRecallsByTf[i], gsPrecisionsByTf[i], recInterpPts)
        macroPrecisions .+= interpPrec
        # ROC Curve   # This is needed for ROC curves. Here Tprs != Recall but interpolated Recall on Fprs
        interpTpr = dedupNInterpolate(gsFprsByTf[i], gsRecallsByTf[i], fprInterpPts)
        macroTprs .+= interpTpr
    end
    
    macroPrecisions ./= length(validTFs)
    macroTprs ./= length(validTFs)
    macroAUPR_interpolated  = sum(0.5 .* (macroPrecisions[2:end] + macroPrecisions[1:end-1]) .* diff(recInterpPts))
    macroAUROC_interpolated = sum(0.5 .* (macroTprs[2:end] + macroTprs[1:end-1]) .* diff(fprInterpPts))

    macroResults = OrderedDict(
        :macroPR => OrderedDict(
            # :auprScalar       => macroAUPR_scalar,
            :auprInterpolated => macroAUPR_interpolated,
            :precisions       => macroPrecisions,
            :recalls          => collect(recInterpPts)
        ),
        :macroROC => OrderedDict(
            # :aurocScalar => macroAUROC_scalar,
            :aurocInterpolated => macroAUROC_interpolated,
            :fprs        => collect(fprInterpPts),
            :tprs        => macroTprs
        )
)
    return macroResults
end

"""
    evaluatePerTF(regs, targs, rankings, infEdges, gsEdgesByTf, uniGsRegs, gsRandPRbyTf;
                  saveDir=nothing, breakTies=true, target_points=1000, min_step=1e-4,
                  step_method=:min_gap, xLimitRecall=1.0)

Compute per-TF precision-recall and ROC metrics for an inferred gene regulatory network (GRN)
against a gold standard, and optionally plot per-TF PR/ROC curves.

# Arguments
- `regs::Vector{<:AbstractString}`       : Regulators for each inferred edge.
- `targs::Vector{<:AbstractString}`      : Target genes for each inferred edge.
- `rankings::Vector{Float64}`            : Confidence scores for inferred edges.
- `infEdges::Vector{String}`             : Inferred edges as `"TF,Target"` strings.
- `gsEdgesByTf::Vector{Array{String}}`   : Gold standard edges grouped by TF.
- `uniGsRegs::Vector{<:AbstractString}`  : Unique TFs in the gold standard.
- `gsRandPRbyTf::Vector{Float64}`      : Random baseline AUPR per TF.

# Keyword Arguments
- `saveDir::Union{String,Nothing}=nothing`  : Directory to save per-TF plots. If `nothing`, plots are skipped.
- `breakTies::Bool=true`                    : Average indicators over tied rankings to smooth PR curves.
- `target_points::Int=1000`                 : Interpolation points for macro curves.
- `min_step::Float64=1e-4`                  : Minimum step size for interpolation.
- `step_method::Symbol=:min_gap`            : Step selection method for macro metrics.
- `xLimitRecall::Float64=1.0`               : : Maximum recall shown on the x-axis of diagnostic per-TF PR plots.

# Returns
An `OrderedDict` with:
- `:perTF`    → `OrderedDict` with `:gsRegs`, `:tfEdges`, `:tfEdgesVec`, `:precisions`,
                `:recalls`, `:fprs`, `:randPR`, `:auprs`, `:arocs`
- `:macroPR`  → `OrderedDict` with `:auprInterpolated`, `:precisions`, `:recalls`
- `:macroROC` → `OrderedDict` with `:aurocInterpolated`, `:fprs`, `:tprs`

# Example
```julia
results = evaluatePerTF(regs, targs, rankings, infEdges, totTargGenes,
                        gsEdgesByTf, uniGsRegs, gsRandPRbyTf)
results[:perTF][:auprs]              # per-TF AUPR values
results[:macroPR][:auprInterpolated] # macro AUPR
```
"""

function evaluatePerTF(
                    regs::Vector{<:AbstractString}, targs::Vector{<:AbstractString}, rankings::Vector{Float64}, infEdges::Vector{String},  
                    gsEdgesByTf::Vector{Array{String}}, uniGsRegs::Vector{<:AbstractString}, gsRandPRbyTf::Vector{Float64};
                    saveDir::Union{AbstractString, Nothing} = nothing, breakTies::Bool=true, target_points=1000, min_step=1e-4, step_method=:min_gap, xLimitRecall::Float64 = 1.0)

    
    if saveDir !== nothing && !isempty(saveDir)
        saveDir = joinpath(saveDir, "perTF")
        mkpath(saveDir)
    end

    totGsRegs = length(uniGsRegs)
    totTargGenes = maximum(length.(gsEdgesByTf))  # approximate max targets

    @info "Computing per-TF metrics"
    gsAuprsByTf = zeros(totGsRegs,1)
    gsArocsByTf = zeros(totGsRegs,1)
    gsPrecisionsByTf = Array{Float64}[]
    gsRecallsByTf = Array{Float64}[]
    gsFprsByTf = Array{Float64}[]  
    
    allTfEdgesVec = Array{Float64}[]  # one vector per TF
    allTfEdges    = Vector{Vector{String}}()  # raw edges
    allTfRanks = Array{Float64}[] # one array per TF
    for (iTF, tf) in enumerate(uniGsRegs)
        # Filter inferred edges for this TF
        inds = findall(x -> x == tf, regs)
        tfEdges = infEdges[inds]
        tfRanks = rankings[inds]
        tfEdgesVec = [in(edge, gsEdgesByTf[iTF]) ? 1 : 0 for edge in tfEdges]
        # Tie-break if needed
        if breakTies
            tfAbsRanks = abs.(tfRanks)
            uniqueRanks = unique(tfAbsRanks)
            meanIndicator = Dict{Float64, Float64}()
            for r in uniqueRanks
                idxs = findall(x -> x == r, tfAbsRanks)
                meanIndicator[r] = mean(tfEdgesVec[idxs])
            end
            tfEdgesVec = [meanIndicator[r] for r in tfAbsRanks]
        end
        # Save for later
        push!(allTfEdgesVec, tfEdgesVec)
        push!(allTfEdges, tfEdges)
        push!(allTfRanks, tfRanks)  # store adjusted ranks for this TF

        # Compute cumulative metrics
        tfTP = 0.0
        tfPrec = Float64[]
        tfRec = Float64[]
        tfFpr = Float64[]
        totNegativesTF = totTargGenes - length(gsEdgesByTf[iTF])
        for j in 1:length(tfEdgesVec)
            tfTP += tfEdgesVec[j]
            fp = j - tfTP
            push!(tfPrec, tfTP / j)
            push!(tfRec, tfTP / length(gsEdgesByTf[iTF]))
            push!(tfFpr, fp / totNegativesTF)
        end
        # Prepend starting point
        if !isempty(tfPrec)
            tfRec = vcat(0.0, tfRec)
            tfPrec = vcat(tfPrec[1], tfPrec)
            tfFpr = vcat(0.0, tfFpr)
        end

        # -- Plot perTF metric
        currRandPR = gsRandPRbyTf[iTF]
        plotPRCurve(tfRec, tfPrec, currRandPR, saveDir; xLimitRecall=xLimitRecall, baseName = tf)
        plotROCCurve(tfFpr, tfRec, saveDir; baseName = tf)

        # Save per-TF metrics
        push!(gsPrecisionsByTf, tfPrec)
        push!(gsRecallsByTf, tfRec)
        push!(gsFprsByTf, tfFpr)
        # AUPR / AUROC per TF
        heightsTF = (tfPrec[2:end] + tfPrec[1:end-1]) / 2
        widthsTF = tfRec[2:end] - tfRec[1:end-1]
        gsAuprsByTf[iTF] = sum(heightsTF .* widthsTF)

        widthsRocTF = tfFpr[2:end] - tfFpr[1:end-1]
        heightsRocTF = (tfRec[2:end] + tfRec[1:end-1]) / 2
        gsArocsByTf[iTF] = sum(widthsRocTF .* heightsRocTF)
    end

    # Compute macro metrics

    # maxLen = maximum(length.(gsPrecisionsByTf))
    # macroPrecisions = zeros(maxLen)
    # macroRecalls = zeros(maxLen)
    # for i in 1:maxLen
    #     precVals = Float64[]
    #     recVals = Float64[]
    #     for (p, r) in zip(gsPrecisionsByTf, gsRecallsByTf)
    #         if i <= length(p)
    #             push!(precVals, p[i])
    #             push!(recVals, r[i])
    #         end
    #     end
    #     macroPrecisions[i] = mean(precVals)
    #     macroRecalls[i] = mean(recVals)
    # end
    # macroAUPR = mean(gsAuprsByTf)
    # macroAUROC = mean(gsArocsByTf)
    
    # A global result but requires perTF to compute
    macroResults = computeMacroMetrics(gsPrecisionsByTf, gsRecallsByTf, gsFprsByTf, gsAuprsByTf, gsArocsByTf;
                    target_points=target_points, min_step=min_step, step_method=step_method)

    perTFDict = OrderedDict(
                :gsRegs => uniGsRegs,
                :tfEdges => allTfEdges,         # raw inferred edges per TF
                :tfEdgesVec => allTfEdgesVec,   # 0/1 indicator per TF
                :precisions => gsPrecisionsByTf,
                :recalls    => gsRecallsByTf,
                :fprs       => gsFprsByTf,
                :randPR     => gsRandPRbyTf,
                :auprs      => gsAuprsByTf,
                :arocs      => gsArocsByTf
        )

    # merge!(perTF, macroResults)
    results = OrderedDict(
        :perTF => perTFDict,
        :macroPR => macroResults[:macroPR],
        :macroROC => macroResults[:macroROC]
    )
    return results
end

"""
    computePR(
        gsFile::String, infTrnFile::String;
        gsRegsFile::Union{String, Nothing} = nothing, targGeneFile::Union{String, Nothing} = nothing,  # filtering
        breakTies::Bool = true, partialAUPRlimit::Float64 = 0.1, doPerTF::Bool = true,                 # computation
        xLimitRecall::Float64 = 1.0, saveDir::Union{String, Nothing} = nothing,                        # plotting
        target_points::Int = 1000, min_step::Float64 = 1e-4, step_method::Symbol = :min_gap)           # interpolation

Compute precision-recall (PR) and ROC metrics for an inferred gene regulatory network (GRN)
against a gold standard, including optional per-TF and macro-level metrics.
Results are saved to a `.jld` file and PR/ROC curves are plotted automatically.

# Arguments
- `gsFile::String`      : Gold standard TSV file with columns in order: TF (1), Target (2), Weight (3) — column names can be anything
- `infTrnFile::String`  : Inferred GRN TSV file with columns in order: TF (1), Target (2), Weight (3) — column names can be anything

# Keyword Arguments
- `gsRegsFile::Union{String,Nothing}=nothing`   : File listing regulators to include. Default uses all GS regulators.
- `targGeneFile::Union{String,Nothing}=nothing` : File listing target genes to include. Default uses all GS targets.
# - `rankColTrn::Int=3`                           : Column index for scores in the inferred GRN file.
- `breakTies::Bool=true`                        : Average indicators over tied rankings to smooth PR curves.
- `partialAUPRlimit::Float64=0.1` : Maximum recall cutoff for computing partial AUPR.
- `xLimitRecall::Float64=1.0`     : Maximum recall shown on the x-axis of diagnostic PR plots.
- `doPerTF::Bool=true`                          : Whether to compute per-TF and macro metrics.
- `saveDir::Union{String,Nothing}=nothing`      : Directory for saving plots and results. Defaults to a timestamped folder.
- `target_points::Int=1000`                     : Interpolation points for macro curves.
- `min_step::Float64=1e-4`                      : Minimum step size for interpolation.
- `step_method::Symbol=:min_gap`                : Step selection method for macro metrics.

# Returns
An `OrderedDict` with:
- `:gsRegs`, `:gsTargs`, `:gsEdges`  → Gold standard regulators, targets, and edges.
- `:randPR`                          → Random baseline PR value.
- `:infRegs`, `:infEdges`            → Inferred GRN regulators and edges.
- `:edgesVec`                        → Tie-adjusted or binary indicator vector.
- `:precisions`, `:recalls`, `:fprs` → Overall PR/ROC curve arrays.
- `:auprs`                           → Dict with `:full` AUPR and `:partial` AUPR (up to `partialAUPRlimit`).
- `:arocs`                           → Overall AUROC.
- `:f1scores`                        → F1-score array.
- `:perTF`                           → Per-TF metrics dict (see `evaluatePerTF`). `nothing` if `doPerTF=false`.
- `:macroPR`                         → Macro PR dict. `nothing` if `doPerTF=false`.
- `:macroROC`                        → Macro ROC dict. `nothing` if `doPerTF=false`.
- `:savedFile`                       → Path to saved `.jld` results file.

# Example
```julia
results = computePR("gs.tsv", "inferredGRN.tsv";
                    breakTies=true, partialAUPRlimit=0.1, saveDir="plots")

results[:auprs][:full]               # overall AUPR
results[:auprs][:partial][:value]    # partial AUPR at recall limit
results[:perTF][:auprs]              # per-TF AUPRs
results[:macroPR][:auprInterpolated] # macro AUPR
```
"""
# function computePR(
#         gsFile::String, infTrnFile::String; 
#         gsRegsFile::Union{String, Nothing} = nothing, targGeneFile::Union{String, Nothing} = nothing,
#         # rankColTrn::Int = 3, 
#         breakTies::Bool = true, partialLimitRecall::Float64 = 0.1, doPerTF::Bool = true,
#         saveDir::Union{String, Nothing} = nothing, target_points::Int = 1000, min_step::Float64 = 1e-4, 
#         step_method::Symbol = :min_gap)
function computePR(
        gsFile::String, infTrnFile::String;
        gsRegsFile::Union{String, Nothing} = nothing, targGeneFile::Union{String, Nothing} = nothing,  # filtering
        breakTies::Bool = true, partialAUPRlimit::Float64 = 0.1, doPerTF::Bool = true,                 # computation
        xLimitRecall::Float64 = 1.0, saveDir::Union{String, Nothing} = nothing,                        # plotting
        target_points::Int = 1000, min_step::Float64 = 1e-4, step_method::Symbol = :min_gap)           # interpolation
    
    # Helper function for empty/fallback PR result
    function emptyPRResult(;
        gsRegs::Vector{String}  = String[],
        gsTargs::Vector{String} = String[],
        gsEdges::Vector{String} = String[],
        infRegs::Vector{String} = String[],
        infEdges::Vector{String} = String[])
        return OrderedDict(
            :gsRegs     => gsRegs,
            :gsTargs    => gsTargs,
            :gsEdges    => gsEdges,
            :infRegs    => infRegs,
            :infEdges   => infEdges,
            :precisions => Float64[],
            :recalls    => Float64[],
            :fprs       => Float64[],
            :auprs      => Dict(
                :full    => 0.0,
                :partial => Dict(
                    :value      => 0.0,
                    :recallLimit => partialAUPRlimit  # captured from outer scope
                )
            ),
            :arocs     => 0.0,
            :f1scores  => Float64[],
            :perTF     => nothing,
            :macroPR   => nothing,
            :macroROC  => nothing,
            :savedFile => nothing
        )
    end
    # ----- Part 1. Load and Validate Input Files
    gsData  = CSV.File(gsFile;     delim="\t", header=true)
    trnData = CSV.File(infTrnFile; delim="\t", header=true)

    validateColumnCount(gsData,  gsFile)
    validateColumnCount(trnData, infTrnFile)

    # Reload selecting only first 3 columns by position
    gsData  = CSV.File(gsFile;     delim="\t", select=[1, 2, 3])
    trnData = CSV.File(infTrnFile; delim="\t", select=[1, 2, 3])


    # Initialize defaults BEFORE any conditionals
    potTargGenes    = String[]
    potTargGenesSet = Set{String}()
    gsPotRegs       = String[]
    gsPotRegsSet    = Set{String}()

    # ----- Part 2. Load Optional Filter Files (Potential Targets and Regulators)
    if targGeneFile !== nothing && !isempty(targGeneFile)
        potTargGenes    = readlines(targGeneFile)
        potTargGenesSet = Set(potTargGenes)
        totTargGenes    = length(unique(potTargGenes))
        @info "Target gene file loaded: $(length(potTargGenes)) genes from $targGeneFile"
    else
        @info "No target gene file provided — will use gold standard targets after GS is loaded"
    end

    # Load regulator list if provided
    if gsRegsFile !== nothing && !isempty(gsRegsFile)
        gsPotRegs    = readlines(gsRegsFile)
        gsPotRegsSet = Set(gsPotRegs)
        @info "Regulator file loaded: $(length(gsPotRegs)) regulators from $gsRegsFile"
    else
        @info "No regulator file provided — will use all regulators in gold standard"
    end

    # ----- Part 3. Filter and Process Gold Standard
    
    # Filter gold standard by regulators and target genes if specified: limit to TF-gene interactions considered by the model
    gsFilteredData = filter(row -> row[3] > 0, gsData)
    # if gsRegsFile !== nothing && !isempty(gsRegsFile)
    #     # gsFilteredData = filter(row -> (row.TF in gsPotRegs) && (row.Target in potTargGenes), gsFilteredData)
    #     gsFilteredData = filter(row ->(row[1] in gsPotRegsSet) && (row[2] in potTargGenesSet), gsFilteredData)
    #     size(gsFilteredData)
    # end
    
    # Filter by regulators if provided
    if !isempty(gsPotRegsSet)
        gsFilteredData = filter(row -> row[1] in gsPotRegsSet, gsFilteredData)
    end

    # Filter by targets if provided  
    if !isempty(potTargGenesSet)
        gsFilteredData = filter(row -> row[2] in potTargGenesSet, gsFilteredData)
    end

    # Edge-case check: skip if no gold standard edges
    if isempty(gsFilteredData)
        @warn "No overlapping gold standard edges after filtering. Skipping PR/ROC evaluation."
        return emptyPRResult()
    end

    # Define gold standard edges (each as "TF,Target")
    gsRegs = collect(row[1] for row in gsFilteredData)
    gsTargs = collect(row[2] for row in gsFilteredData)
    totGsInts = size(gsFilteredData)[1]
    uniGsRegs = unique(gsRegs) # unique regulators or TFs in GS
    totGsRegs = length(uniGsRegs)
    uniGsTargs = unique(gsTargs) # unique targets in GS
    totGsTargs = length(uniGsTargs)
    # Create new edges vector. (Each edge is a tuple (TF, gene))
    gsEdges = [string(gsRegs[i], ",", gsTargs[i]) for i in 1:length(gsRegs)]

    # If targGeneFile wasn’t provided, use GS targets.
    if targGeneFile === nothing || isempty(targGeneFile)
        potTargGenes = uniGsTargs
        totTargGenes = length(potTargGenes)
    end

    # Compute evaluation universe and random PR baseline
    gsTotPotInts = totTargGenes*totGsRegs  # complete universe size
    # gsTotPotInts = totGsRegs * totGsTargs
    gsRandPR = totGsInts/gsTotPotInts
    
    # group gold standard edges by regulator
    gsEdgesByTf = Array{String}[]
    gsRandPRbyTf = zeros(totGsRegs)
    for gind = 1:totGsRegs
        currInds = findall(x -> x == uniGsRegs[gind], gsRegs)
        push!(gsEdgesByTf, vec(permutedims(gsEdges[currInds])))
        gsRandPRbyTf[gind] = length(gsEdgesByTf[gind])/totGsTargs
    end

    # ----- Part 4. Load and Filter Inferred GRN
    @info "Loading and Processing Inferred GRN"
    # Filter "inferred GRN" to include only 'TFs' in GS and 'TargetGenes' in 'potTargGenes'
    grnData = filter(row -> (row[1] in uniGsRegs) && (row[2] in potTargGenes), trnData)
    grnData = collect(grnData)

    if isempty(grnData)
        @warn "No inferred edges after filtering. Skipping evaluation."
        return emptyPRResult(gsRegs=uniGsRegs, gsTargs=uniGsTargs, gsEdges=gsEdges)
    end

    # Order by absolute value of weights/confidences
    grnData = sort(grnData, by =  row -> abs(row[3]), rev = true)
    regs = collect(row[1] for row in grnData)
    targs = collect(row[2] for row in grnData)
    rankings = collect(row[3] for row in grnData)
    absRankings = abs.(rankings)
    # Create inferred edges vector. (Each edge is a tuple (TF, gene))
    infEdges = [string(regs[i], ",", targs[i]) for i in 1:length(regs)] 
    totTrnInts = length(infEdges)

    # ----- Part 5. Compute Edge Indicators
    @info "Computing Edge Indicators and Labels"
    # create binary labels. 1 if infEgde in gsEdge and 0 otherwise
    commonEdgesBinaryVec = [in(edge, gsEdges) ? 1 : 0 for edge in infEdges]
    @info "Total Interactions w/ GS TFs ($(length(unique(regs)))): $totTrnInts"

    if breakTies   # If breaking ties is desired, compute tie-adjusted (mean) vector.
        # Break ties in weights/confidences to smooth out abrupt jumps caused by abitrary ordering of tied predictions
        # Create a dictionary mapping each score to the mean value of commonEdgesVec for tied predictions
        @info "Tie breaker enabled"
        uniqueRankings = unique(absRankings)
        meanIndicator = Dict{Float64, Float64}()
        for currRank in uniqueRankings
            inds = findall(x -> x == currRank, absRankings)
            meanIndicator[currRank] = mean(commonEdgesBinaryVec[inds])
        end
        meanEdgesVec = [meanIndicator[ix] for ix in absRankings] 
        edgesVec = meanEdgesVec
    else  # If not breaking ties, use binary vector
        edgesVec = commonEdgesBinaryVec
    end
    
    # ----- Part 6. Compute Overall PR/ROC Performance Metrics
    @info "Computing Performance Metrics"
    totalNegatives = gsTotPotInts - totGsInts   # Total posisble interactions - length(gsEdges) = TN + FP
    gsPrecisions = zeros(totTrnInts)
    gsRecalls = zeros(totTrnInts)
    gsFprs = zeros(totTrnInts)

    cummulativeTP = 0.0   # can be fractional in tie-adjusted mode
    for idx in 1:totTrnInts
        cummulativeTP +=  edgesVec[idx]  # Add the effective contirbution for this prediction. This is also True Positive
        # False positive (FP). This is a weighted FP in the case of tie breaking
        falsePositives = idx - cummulativeTP

        # Compute precision: TP divided by total prediction so far
        # idx is always the total predicted positive at any point j 
        # such that j = TP + FP
        gsPrecisions[idx] = cummulativeTP / idx # 
        gsRecalls[idx] = cummulativeTP / totGsInts  # truePositives/length(gsEdges) == tp/tp+fn
        gsFprs[idx] = falsePositives / totalNegatives
    end
    
    # Prepend starting point for plotting.
    if !isempty(gsPrecisions)
        gsRecalls = vcat(0.0, gsRecalls)
        gsPrecisions = vcat(gsPrecisions[1], gsPrecisions)
        gsFprs = vcat(0.0, gsFprs)
    end
    # Compute F1-scores
    # gsF1scores = 2 * (gsPrecisions .* gsRecalls) ./ (gsPrecisions + gsRecalls);
    gsF1scores = [p + r > 0 ? 2p*r/(p+r) : 0.0 for (p, r) in zip(gsPrecisions, gsRecalls)];

    # ----- Part 7. Compute AUPR and AUROC
    @info "Computing AUPR and AUROC"
    # Here, AUPR is computed using trapezoidal rule. Other methods available is a step-function approximation
    heights = (gsPrecisions[2:end] + gsPrecisions[1:end-1])/2  
    widths = gsRecalls[2:end] - gsRecalls[1:end-1]
    gsAuprs = sum(heights .* widths)
    #= 
    # Step-function approximation. This works but is less robust
    gsAuprs = 0.0
    prev_recall = 0.0
    for (r, p) in zip(gsRecalls, gsPrecisions)
        delta = r - prev_recall
        gsAuprs += delta * p
        prev_recall = r
    end
    =#

    # PARTIAL AUPR @ partialAUPRlimit
    indx = findall(r -> r <= partialAUPRlimit, gsRecalls) # Find indices of recalls <= partialUpperLimRecall
    recalls_sub = gsRecalls[indx]
    precisions_sub = gsPrecisions[indx]
    heights_sub = (precisions_sub[2:end] + precisions_sub[1:end-1]) ./ 2
    widths_sub = recalls_sub[2:end] - recalls_sub[1:end-1]
    partialAUPR = sum(heights_sub .* widths_sub)

    # AROC : Trapezoidal Rule
    widthsRoc = gsFprs[2:end] - gsFprs[1:end-1]    # Change in FPR (the x-axis) between successive points.
    heightsRoc = (gsRecalls[2:end] + gsRecalls[1:end-1]) / 2  # Average TPR (recall) for each segment.
    gsArocs = sum(widthsRoc .* heightsRoc)
    #=
    gsArocs = 0.0
    prev_fpr = 0.0
    for (f, r) in zip(fprs, recalls)
        delta = f - prev_fpr
        gsArocs += delta * r
        prev_fpr = f
    end
    =#

    # ----- Part 8. Save Directory
    baseName = splitext(basename(infTrnFile))[1]
    if saveDir === nothing || isempty(saveDir)
        dateStr  = Dates.format(now(), "yyyymmdd_HHMMSS")
        saveDir  = joinpath(pwd(), baseName * "_" * dateStr)
    end
    mkpath(saveDir)

    # ----- Part 9. Per-TF and Macro Metrics
    resultsTF   = nothing
    perTFDict   = nothing
    macroPRDict = nothing
    macroROCDict = nothing
    if doPerTF
        if !isempty(regs) && !isempty(targs) && !isempty(gsEdgesByTf)
            resultsTF = evaluatePerTF(
                regs, targs, rankings, infEdges, gsEdgesByTf, uniGsRegs, gsRandPRbyTf;
                saveDir,
                breakTies=breakTies, target_points=target_points, min_step=min_step,
                step_method=step_method, xLimitRecall=xLimitRecall
            )
            # Extract Per-TF Results
            perTFDict    = resultsTF[:perTF]
            macroPRDict  = resultsTF[:macroPR]
            macroROCDict = resultsTF[:macroROC]
        else
            @warn "Skipping per-TF evaluation: no valid edges after filtering"
        end
    end


    # ----- Part 10. Plot and Save Results
    suffix = breakTies ? "_tiesBroken" : ""
    baseName = baseName * suffix
    plotPRCurve(gsRecalls, gsPrecisions, gsRandPR, saveDir; xLimitRecall=xLimitRecall, baseName = baseName)
    plotROCCurve(gsFprs, gsRecalls, saveDir; baseName = baseName)


    results = OrderedDict(
        :gsRegs => uniGsRegs,
        :gsTargs => uniGsTargs,
        :gsEdges => gsEdges,
        :randPR => gsRandPR,
        :infRegs => regs,
        :infEdges => infEdges,
        :breakTies => breakTies,
        :stepVals => rankings,
        :rankings => rankings,
        # :commonEdgesBinaryVec => commonEdgesBinaryVec,
        # :meanEdgesVec => breakTies ? meanEdgesVec : nothing,
        :edgesVec => edgesVec,
        :precisions => gsPrecisions,
        :recalls => gsRecalls,
        :fprs => gsFprs,
        # :auprs => gsAuprs,
        :auprs => Dict(
            :full => gsAuprs,
            :partial => Dict(
                :value => partialAUPR,
                :recallLimit => partialAUPRlimit
            )
        ),
        :arocs => gsArocs,
        :f1scores => gsF1scores,
        # --- Include results from evaluatePerTF ---
        :perTF => perTFDict,
        :macroPR => macroPRDict,
        :macroROC => macroROCDict   
    )
    
    # Define a standard filename for saving the performance metrics
    savedFile = joinpath(saveDir, baseName * "_PerformanceMetric.jld")
    @save savedFile results

    # Add the saved file path to the results, so the caller immediately knows it
    results[:savedFile] = savedFile

    return results
end