"""
    GRN

This module provides functionality for building and combining Gene Regulatory Networks (GRNs).

The GRN files are assumed to have the following columns:
                        TF, Gene, signedQuantile, Stability, Correlation, strokeVals, strokeWidth, inPrior
Main function:
- `combineGRNs`: combine multiple GRN files with options for mean/max aggregation, 
  controlling edges per gene, and saving the combined network.

Dependencies:
- Uses `PriorUtils` for prior-related helper functions.
- Includes internal utilities in `utilsGRN.jl`.
"""

# using ..DataUtils
# using CSV, DataFrames, Statistics, StatsBase, Printf, Dates
# using ArgParse

# export combineGRNs
# ────────────────────────────────────────────────────────────────
#  Helper: deterministic tie-breaker
# ────────────────────────────────────────────────────────────────
function pickRow(g, primary::Union{String,Symbol}; order::Symbol = :max)

    """
    pickRow(g, primary; order = :max)

    Return a single row index from SubDataFrame `g` according to a deterministic
    tie-breaking hierarchy.

    Arguments
    ─────────
    g       - SubDataFrame produced by `groupby`
    primary - "Stability"/:Stability  or  "Quantile"/:Quantile
    order   - :max  (default) or :min   (optimise up or down)

    Tie-breaking
    ────────────
    If primary = Stability   :  Stability → Quantile → |Correlation|
    If primary = Quantile    :  Quantile  → |Correlation|
    """

    prim = Symbol(primary)                  # normalise to Symbol

    # 1) extreme value of the primary column
    pVals   = g[!, prim]
    extreme = order === :max ? maximum(pVals) : minimum(pVals)
    idxs    = findall(==(extreme), pVals)
    length(idxs) == 1 && return idxs[1]

    # 2) additional tie-breakers
    if prim == :Stability && "signedQuantile" in names(g)
        qVals = g.signedQuantile[idxs]
        extQ  = order === :max ? maximum(qVals) : minimum(qVals)
        idxs  = idxs[findall(==(extQ), qVals)]
        length(idxs) == 1 && return idxs[1]
        # if still tied we fall through to correlation
    end

    if "Correlation" in names(g)
        absCorr = abs.(g.Correlation[idxs])
        bestC   = maximum(absCorr)                      # always take largest |ρ|
        idxs    = idxs[findall(==(bestC), absCorr)]
    end

    return idxs[1]      # deterministic fallback
end


# ────────────────────────────────────────────────────────────────
# Main function
# ────────────────────────────────────────────────────────────────
function aggregateNetworks(nets2combine::Vector{String}; combineOpt::Union{String,Symbol}, meanEdgesPerGene::Int,  useMeanEdgesPerGeneMode::Bool, 
                    saveDir::Union{String,Nothing}=nothing, saveName::String="")

    """
    # Arguments
    - `nets2combine::Vector{String}`: A vector of GRN file names to combine.
    - `combineOpt::Union{String,Symbol}`: Specifies the combination strategy (e.g., "max", "meanQuantile").
    - `meanEdgesPerGene::Int`: The mean number of edges per gene.
    - `saveDir::Union{String,Nothing}`: Optional. Directory to save results. Defaults to `nothing` (do not save).
    - `saveName::String`: Optional. Filename for saving combined GRN. Defaults to `""`.
        useMeanEdgesPerGene = false # This option controls whether edge selection is done per-group or globally. If `true`, selects the top `meanEdgesPerGene` edges **per target gene**. 
                                # If `false`, selects a flat total number of edges equal to `length(unique(targs)) * meanEdgesPerGene`.
    """
    allDfs = DataFrame[]
    
    # ──── 1. Read each file and then combine; create a rank column if not present ───────────────────────────
    for file in nets2combine
        # Adjust CSV reading as needed (TSV expected if tab sep)
        df = CSV.read(file, DataFrame; delim='\t')
        # println(size(df))
        
        # # Ensure the required columns exist
        # requiredCols = ["TF", "Gene", "signedQuantile", "Stability", "Correlation", "inPrior"]
        # dfCols = names(df)
        # for col in requiredCols
        #     if !(col in dfCols)
        #         error("Column $(col) not found in file $(file)")
        #     end
        # end 
        # df = df[!, requiredCols]

        # ensure numeric columns are Float64
        for col in [:Stability, :Correlation]
            if !(eltype(df[!, col]) <: AbstractFloat)
                df[!, col] = parse.(Float64, string.(df[!, col]))
            end
        end
        
        # fill in ranks ONLY if these columns do not already exist
        if !("signedQuantile" in names(df)) #|| !(:Ranks in names(df))
            st_ecdf = ecdf(df.Stability)
            quantile = st_ecdf.(df.Stability)
            df.signedQuantiles  =  quantile .* sign.(df.Correlation)
        end
        push!(allDfs, df)
    end

    # ──── 2. CONCATENATE AND AGGREGATE ACCORDING TO combineOpt ───────────────────────────
    # Combine all dataframes vertically
    combinedDf = vcat(allDfs...)
    # Group by TF and Gene and aggregate according to the chosen method.
    groupedDf = groupby(combinedDf, ["TF", "Gene"])
    aggRows = NamedTuple[]         

    for g in groupedDf
        # All rows in this group share the same TF, Gene and inPrior
        tf, gene, dash = g.TF[1], g.Gene[1], g.inPrior[1]
        # tf, gene, dash = g.TF[1], g.Gene[1], g.strokeDashArray[1]

        if combineOpt == "mean"
            c = mean(g.Correlation)
            push!(aggRows, (TF = tf, Gene = gene, Stability =  mean(g.Stability), 
                            Correlation = c,
                            inPrior = dash))

        elseif combineOpt == "max"
            # idx = combineOpt == "max" ? argmax(g.Stability) : argmin(g.Stability)
            idx = pickRow(g, "Stability", order = :max)
            row = g[idx, :]
            push!(aggRows, (TF = tf, Gene = gene,  signedQuantile =  row.signedQuantile, Stability = row.Stability,
                        Correlation = row.Correlation, inPrior = row.inPrior))
        
        elseif combineOpt == "min"
            # idx = combineOpt == "max" ? argmax(g.Stability) : argmin(g.Stability)
            idx = pickRow(g, "Stability", order = :min)
            row = g[idx, :]
            push!(aggRows, (TF = tf, Gene = gene,  Stability = row.Stability,
                            Correlation = row.Correlation, inPrior = row.inPrior))

        elseif combineOpt == "meanQuantile"
            c = mean(g.Correlation)
            push!(aggRows, (TF = tf, Gene = gene, signedQuantile =  mean(g.signedQuantile), Stability = mean(g.Stability),
                            Correlation = c,
                            inPrior = dash))

        elseif combineOpt == "maxQuantile"   
            # idx = argmax(g.signedQuantile)
            idx = pickRow(g, "signedQuantile", order = :max)
            row = g[idx, :]
            push!(aggRows, (TF = tf, Gene = gene, signedQuantile = row.signedQuantile, Stability = row.Stability, Correlation = row.Correlation,
                            inPrior = row.inPrior))
        else
            error("Invalid combineOpt: $(combineOpt).")
        end
    end

    aggregatedDf = DataFrame(aggRows)
    

    # ────── 3. Filter out extrememly weak correlations and compute new signed quantiles using broadcasting ──────────────────────
    pcut = 0.01 # Correlation cutoff
    aggregatedDf = filter(:Correlation => x -> abs(x) > pcut, aggregatedDf)
    # ── 3b. Compute signedQuantile  ─────────────────
    if combineOpt in ("max","min","mean")
        stability = aggregatedDf.Stability          # Vector
        F = ecdf(stability)                         # F(x) = proportion ≤ x
        signed_q = F.(stability) .* sign.(aggregatedDf.Correlation)
        aggregatedDf[!,:signedQuantile] = signed_q
    end
    
    # Checks:
    # abs_q = abs.(aggregatedDf.signedQuantile)
    # i_min = argmin(abs_q)
    # row_min = aggregatedDf[i_min, :]   

    # sort for reproducibility
    sort!(aggregatedDf, :signedQuantile, by= abs, rev  = true)
    # sort!(aggregatedDf, combineOpt in ("meanQuantile", "maxQuantile") ? :signedQuantile : :Stability,
    #     by  = combineOpt in ("meanQuantile", "maxQuantile") ? abs : identity, rev = true)

    # ────── 5. choose top meanEdgesPerGene * nGenes edges ──────────────────────

    if useMeanEdgesPerGeneMode
        println("Selecting the top (meanEdgesPerGene * unique targets) edges")
        # Select the top `topEdgeCount = meanEdgesPerGene * uniqueGenes` & compute quantiles for all rankings
        uniqueGenes = length(unique(aggregatedDf.Gene))   # current number of unique targets
        topEdgeCount = round(Int, meanEdgesPerGene * uniqueGenes)     # Compute the number of edges to select
        selectionIndices = 1:min(topEdgeCount, nrow(aggregatedDf))
    else
        println("Selecting the top meanEdgesPerGene edges per target gene")
        selectionIndices = firstNByGroup(aggregatedDf.Gene, meanEdgesPerGene)
    end
    selectedDf = aggregatedDf[selectionIndices, :]
    println("Selected $(length(selectionIndices)) edges using meanEdgesPerGene = $meanEdgesPerGene.")

    # # ────── 6. Generate strokeWidth and colors
    # # ────── Color calculations for JP-Gene-Viz 
    # minRank, maxRank = extrema(selectedDf.Stability)
    # rankRange = max(maxRank - minRank, eps())  # prevent division by zero
    # strokeWidth = 1 .+ (selectedDf.Stability .- minRank) ./ rankRange
    # # Color mapping
    # medRed = [228, 26, 28]
    # lightGrey = [217, 217, 217]
    # strokeVals = map(abs.(selectedDf.Correlation)) do corr
    #         color = corr * medRed .+ (1 - corr) * lightGrey
    #         "rgb($(floor(Int, round(color[1]))),$(floor(Int, round(color[2]))),$(floor(Int, round(color[3]))))"
    #         # "rgb(" * string(floor(Int, round(color[1]))) * "," * string(floor(Int, round(color[2]))) * "," * string(floor(Int, round(color[3]))) * ")"
    # end
    # selectedDf[!, :strokeVals] = strokeVals
    # selectedDf[!, :strokeWidth] = strokeWidth


    # ── 7. save to disk if saveDir is provided ──────────────────────────────────────────
    if saveDir !== nothing
        mkpath(saveDir)
        tag = isnothing(saveName) || isempty(saveName) ? combineOpt : "$(saveName)_$(combineOpt)"
        # CSV.write(joinpath(saveDir, "$(tag)_aggregated.tsv"), aggregatedDf; delim = '\t')
        CSV.write(joinpath(saveDir, "combined_$(tag).tsv"),   selectedDf;  delim = '\t')

        wideDf = convertToWide(selectedDf[ :,["TF","Gene","signedQuantile"]]; indices=(2,1,3))
        wideDf = coalesce.(wideDf, 0)
        wideFile = joinpath(saveDir, "combined_$(tag)_sp.tsv")
        writeTSVWithEmptyFirstHeader(wideDf, wideFile; delim = '\t')
        println("Files written under ", saveDir)
    end

    return selectedDf
end

# -----------------------
# main procedure
# -----------------------
# This section executes if the script is run directly. It will not execute if called/imported into another script.
if abspath(PROGRAM_FILE) == @__FILE__
    using ArgParse

    # Define argument parser
    s = ArgParseSettings()
    @add_arg_table s begin
        "--combineOpt"
            help = "Combination option: max, min, or mean (default: mean)"
            arg_type = String
            default = "mean"
        "--meanEdgesPerGene"
            help = "Mean number of edges per unique gene."
            arg_type = Int
            required = true
        "--useMeanEdgesPerGeneMode"
            help = "controls whether edge selection is done per-group or globally. If true, selects length(unique(targs)) * meanEdgesPerGene edges.
                    If `false`, selects selects the top meanEdgesPerGene edges per target gene." 
            arg_type = Bool
            default = true
        "--saveDir"
            help = "Directory in which to save output files."
            arg_type = String
            default = ""
        "--saveName"
            help = "Base name for the saved file."
            arg_type = String
            default = ""
        "files..."
            help = "List of GRN files (TSV format) to combine."
    end

    parsedArgs = parse_args(s)
    fileList = parsedArgs["files"]
    combineOpt = parsedArgs["combineOpt"]
    meanEdgesPerGene = parsedArgs["meanEdgesPerGene"]
    useMeanEdgesPerGeneMode = parsedArgs["useMeanEdgesPerGeneMode"]
    saveDir = isempty(parsedArgs["saveDir"]) ? nothing : parsedArgs["saveDir"]
    saveName = parsedArgs["saveName"]

    @printf("Combining %d files with combination option: %s\n", length(fileList), combineOpt)

    selectedDf = aggregateNetworks(fileList; combineOpt=combineOpt, meanEdgesPerGene=meanEdgesPerGene, saveDir=saveDir, saveName=saveName)

    # Optionally, print summary
    @printf("Final network has %d edges spanning %d unique genes.\n", nrow(selectedDf), length(unique(selectedDf.Gene)))
end

