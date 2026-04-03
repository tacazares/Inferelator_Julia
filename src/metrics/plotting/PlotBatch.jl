# ----------------------------------------------------------------------
# Color Handling
# ----------------------------------------------------------------------

"""
    padColors(listFiles)

Ensures enough colors for plotting by extending a predefined palette.
Returns a vector of hex color codes, generating random colors if needed.

# Arguments
- `listFiles`: Vector of items needing unique colors.
"""
function padColors(listFiles)
    lineColors = copy(DEFAULT_COLORS) 
    lenFile = length(listFiles)
    lenColor =  length(lineColors)
    excessCT = lenFile - lenColor

    if excessCT > 0
        @warn "More entries than available colors — generating random colors for $excessCT excess entries"
        colorGen = [hex(RGB(rand(), rand(), rand())) for _ in 1:(excessCT + 12)]  # Generate random colors for plots
        uniqueCols = setdiff(colorGen, lineColors)
        # Pad lineColors with the unique additional colors
        lineColors = vcat(lineColors, uniqueCols[1:excessCT])
    end
    return lineColors
end 

# ----------------------------------------------------------------------
# Helper: Data Loading
# ----------------------------------------------------------------------
# function loadPRData(source)
#     if isa(source, String)
#         try
#             data = load(source)
#             return Dict(
#                 :precisions => data["results"][:precisions],
#                 :recalls    => data["results"][:recalls],
#                 :randPR     => get(data["results"], :randPR, nothing)
#             )
#         catch e
#             @warn "Could not load $source" exception=(e, catch_backtrace())
#             return nothing
#         end
#     elseif isa(source, Dict) || isa(source, OrderedDict)
#         if haskey(source, :precisions) && haskey(source, :recalls)
#             return Dict(:precisions=>source[:precisions], :recalls=>source[:recalls],
#                         :randPR=>get(source, :randPR, nothing))
#         elseif haskey(source, "results")
#             r = source["results"]
#             return Dict(:precisions=>r[:precisions], :recalls=>r[:recalls],
#                         :randPR=>get(r, :randPR, nothing))
#         end
#     end
#     return nothing
# end

function loadPRData(source; mode::Symbol=:global)
    # Load file safely if source is a path
    r = if isa(source, String)
        try
            load(source)["results"]
        catch e
            @warn "Could not load $source" exception=(e, catch_backtrace())
            return nothing
        end
    else
        source
    end

    if r === nothing
        return nothing
    end

    # Determine mode
    if mode == :macro
        if !haskey(r, :macroPR)
            @warn "No macroPR found, falling back to global"
            r_macro = r
        else
            r_macro = r[:macroPR]
        end
        return Dict(
            :precisions => r_macro[:precisions],
            :recalls    => r_macro[:recalls],
            :randPR     => get(r_macro, :randPR, get(r, :randPR, nothing))  # fallback to global randPR
        )
    elseif mode == :global
        return Dict(
            :precisions => r[:precisions],
            :recalls    => r[:recalls],
            :randPR     => get(r, :randPR, nothing)
        )
    elseif mode == :perTF
        if !haskey(r, :perTF)
            @warn "No perTF data found"
            return nothing
        end
        return r[:perTF]
    else
        error("Unknown mode: $mode")
    end
end

# ----------------------------------------------------------------------
# Transformation Helper
# ----------------------------------------------------------------------
"""
    function makeTransform(yScale::String)
Creates forward and inverse transformation functions for plot scaling.

# Arguments
- `yScale::String`: Transformation type ("sqrt", "cubert", or "linear")

Returns a tuple of (forward, inverse) functions. Defaults to identity functions if scale type not recognized.
    """
function makeTransform(yScale::String)
    if yScale == "linear"
        forwardFunc = x -> x
        inverseFunc = x -> x
    elseif yScale == "sqrt"
        forwardFunc = x -> sqrt.(x)
        inverseFunc = x -> x .^ 2
    elseif yScale == "cubert"
        forwardFunc = x -> x .^ (1/3)
        inverseFunc = x -> x .^ 3
    else
        # For "linear" or any unknown yScale, return identity functions
        @warn "Unknown scale type: $yScale. Using linear scale."
        forwardFunc = x -> x
        inverseFunc = x -> x
    end
    return forwardFunc, inverseFunc
end

# ----------------------------------------------------------------------
# Plotting Helpers
# ----------------------------------------------------------------------

# # Helper function: Plot data on one or more axes.
# function plotFileData!(axes, legendLabel::String, dataSource, color; lineType=nothing)
#     """
#     Plots precision-recall curves on specified axes from data file or direct data.

#     # Arguments
#     - `axes`: Array of plot axes to draw on
#     - `legendLabel::String`: Label for plot legend
#     - `dataSource`: Either a file path (String) or a Dict containing PR data
#     - `color`: Color for plotting
#     - `lineType=nothing`: Optional line style specification

#     # Returns
#     Loaded/processed data or nothing if processing fails.
#     """
    
#     # Initialize variables to hold the data
#     precisions = nothing
#     recalls = nothing
#     randPR = nothing
#     fileData = nothing
    
#     # Process the data source based on its type
#     if isa(dataSource, String)
#         # It's a file path - try to load it
#         try
#             fileData = load(dataSource)
#             # Extract data from the loaded file
#             precisions = fileData["results"][:precisions]
#             recalls = fileData["results"][:recalls]
#             randPR = get(fileData["results"], :randPR, nothing)
#         catch e
#             println("Error loading file: $dataSource - $e")
#             return nothing
#         end
#     elseif isa(dataSource, Dict) || isa(dataSource, OrderedDict)
#         # It's already a data dictionary
#         fileData = dataSource
        
#         # Check if it's a results dictionary or a direct data dictionary
#         if haskey(dataSource, :precisions) && haskey(dataSource, :recalls)
#             # Direct data format
#             precisions = dataSource[:precisions]
#             recalls = dataSource[:recalls]
#             randPR = get(dataSource, :randPR, nothing)
#         elseif haskey(dataSource, "results")
#             # Nested results format (like from a loaded JLD file)
#             precisions = dataSource["results"][:precisions]
#             recalls = dataSource["results"][:recalls]
#             randPR = get(dataSource["results"], :randPR, nothing)
#         else
#             println("Error: Data dictionary does not contain required precision/recall data")
#             return nothing
#         end
#     else
#         println("Error: Unsupported data source type: $(typeof(dataSource))")
#         return nothing
#     end
    
#     # Ensure we have valid data before plotting
#     if isnothing(precisions) || isnothing(recalls)
#         println("Error: Could not extract precision/recall data from source")
#         return nothing
#     end
    
#     # Plot the data on each provided axis
#     for ax in axes
#         if isnothing(lineType) || lineType == ""
#             # Use default linestyle
#             ax.plot(recalls, precisions, label=legendLabel, color=color, linewidth=0.8)
#         else
#             # Specify the provided line type (linestyle)
#             ax.plot(recalls, precisions, label=legendLabel, color=color, linewidth=0.5, linestyle=lineType)
#         end
#     end
    
#     return Dict(
#         # "data" => fileData,
#         "precisions" => precisions,
#         "recalls" => recalls,
#         "randPR" => randPR
#     )
# end

function plotFileData!(axes::AbstractVector, source, label::String, color::String; lineType=nothing, lineWidth=nothing, mode::Symbol=:global)
    data = loadPRData(source; mode)
    if isnothing(data)
        @warn "Skipping $label ($source) — no valid PR data found"
        return false
    end

    for ax in axes
        ax.plot(data[:recalls], data[:precisions], label=label, color=color,linewidth=(lineWidth === nothing ? 0.7 : lineWidth),
        linestyle=(lineType === nothing ? "-" : lineType))
    end
    return data
end


# Function to combine legend handles and labels from every axis.
function combineLegends(fig)
    allHandles = Any[]
    allLabels = String[]
    # Loop over each axis in the figure.
    for ax in fig.axes
        handles, labels = ax.get_legend_handles_labels()
        for (h, l) in zip(handles, labels)
            if !(l in allLabels)
                push!(allHandles, h)
                push!(allLabels, l)
             end
        end
    end
    return allHandles, allLabels
end

function styleAxis!(ax; xZoom=0.1, yRange::Union{Nothing,Tuple}=nothing,
    xStep=0.05, yStep=0.1, yScale="linear", setXticks=true, tickLabelSize::Int = 7)
# function styleAxis!(ax; xZoom, yRange, xStep, yStep, yScale, setXticks, tickLabelSize)
    # X-axis
    ax.set_xlim(0, xZoom)
    if setXticks
        xMax = mod(xZoom, xStep) == 0 ? xZoom : ceil(xZoom/xStep) * xStep
        ax.set_xticks(0:xStep:xMax)
    end

    # Y-axis
    forwardFunc, inverseFunc = makeTransform(yScale)
    yMin, yMax = isnothing(yRange) ? (0.0, 1.0) : yRange
    ax.set_ylim(yMin, yMax)

    # Y ticks
    yTicks = yMin:yStep:yMax
    ax.set_yticks(yTicks)
    ax.set_yticklabels(string.(round.(yTicks; digits=2)))
    ax.set_yscale("function", functions=(forwardFunc, inverseFunc))

    # Labels
    # ax.set_xlabel(xlabel, fontsize=axisTitleSize, family="Nimbus Sans")
    # ax.set_ylabel(ylabel, fontsize=axisTitleSize, family="Nimbus Sans")

    # Grid & ticks
    ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
    ax.minorticks_on()
    ax.grid(true, which="minor", linestyle=":", linewidth=0.25, color="lightgray")
    ax.tick_params(axis="both", which="both", labelsize=tickLabelSize, direction="out")
end

"""
    plotPRCurves(listFilePR, dirOut, saveName; kwargs...)

Plot and save precision-recall curves for one or more networks against a gold standard.
Supports single-axis and broken y-axis modes. Saves both legend and no-legend versions.

# Arguments
- `listFilePR`      : `OrderedDict` mapping legend labels to either file paths or
                      `OrderedDict`s with `:precisions`, `:recalls`, and `:randPR` keys.
- `dirOut::String`  : Output directory for saved plots.
- `saveName::String`: Base name for output files.

# Keyword Arguments
- `xLimitRecall::Float64=0.1`           : Maximum recall shown on x-axis.
- `yZoomPR::Vector=[]`                  : Y-axis display range.
                                          `[]`         → full range 0–1.
                                          `[0.8]`      → limit y-axis to 0–0.8.
                                          `[0.4, 0.9]` → broken y-axis with gap between 0.4 and 0.9.
- `xStepSize::Union{Nothing,Real}=nothing` : X-axis tick step size. Default auto-set to 0.05.
- `yStepSize::Union{Nothing,Real}=nothing` : Y-axis tick step size. Default auto-set to 0.2.
- `yScale::String="linear"`             : Y-axis scale — `"linear"`, `"sqrt"`, or `"cubert"`.
- `isInside::Bool=true`                 : Place legend inside (`true`) or outside (`false`) the plot.
- `lineColors::Vector=[]`               : Custom line colors. Empty uses default palette.
- `lineTypes::Vector=[]`                : Custom line styles. Empty uses solid lines.
- `lineWidths::Vector=[]`               : Custom line widths. Empty uses default width.
- `heightRatios=nothing`                : Height ratios for broken y-axis panels `[topHeight, bottomHeight]`.
                                        `nothing` --> auto-computed proportionally from `yZoomPR` values.
                                        Pass explicit values e.g. `[3.0, 0.5]` to control panel sizes manually —
                                        larger value = taller panel.
- `mode::Symbol=:global`                : PR data to use:
                                          `:global` → overall PR curve.
                                          `:macro`  → macro-averaged PR (falls back to global if missing).
                                          `:perTF`  → per-TF PR curves.

# Returns
Path to the saved plot file.

# Example
```julia
plotPRCurves(listFilePR, "plots", "MyNetwork";
             xLimitRecall=0.1, yZoomPR=[0.4, 0.9],
             isInside=false, mode=:macro)
```
"""
function plotPRCurves(listFilePR, dirOut::String, saveName::String;
                    xLimitRecall = 0.1, yZoomPR = [], xStepSize::Union{Nothing,Real}=nothing,
                    yStepSize::Union{Nothing,Real}=nothing,
                    yScale::String="linear", isInside::Bool=true,
                    lineColors=[], # empty vector means "use default"
                    lineTypes=[], # empty vector means "use default",
                    lineWidths = [],
                    heightRatios=[0.5, 3.0],
                    mode::Symbol = :global 
                    )

    # Defaults 
    xStep = isnothing(xStepSize) ? 0.05 : xStepSize
    yStep = isnothing(yStepSize) ? 0.2  : yStepSize
    dirOut = isempty(dirOut) ? pwd() : mkpath(dirOut)
    # Style parameters
    axisTitleSize = 9
    tickLabelSize = 7
    # plotTitleSize = 16
    legendSize = 9

    # Color assignment
    if isempty(lineColors)
        lineColors = padColors(listFilePR)
    end

    # Auto height ratios if broken axis
    if length(yZoomPR) == 2 && heightRatios === nothing
        lowerSpan = max(yZoomPR[1] - 0.0, 0.0)      # bottom panel span (0 -> yZoomPR[1])
        upperSpan = max(1.0 - yZoomPR[2], 0.0)      # top panel span (yZoomPR[2] -> 1)
        # protect against zero spans
        if upperSpan == 0.0
            upperSpan = 1e-6
        end
        if lowerSpan == 0.0
            lowerSpan = 1e-6
        end
        # gridspec expects [top_height, bottom_height]
        heightRatios = [upperSpan, lowerSpan]
    end
    # lastPlotData = nothing   # will store last loaded file data (for randPR)
    # println(">>> Inside plotPRCurves: yZoomPR = ", yZoomPR, "  length=", length(yZoomPR))

    ## Making Plots
    if isempty(yZoomPR) || length(yZoomPR) == 1
        yZoom = isempty(yZoomPR) ? 1.0 : yZoomPR[1]

        # Create a Single PyPlot figure
        fig, ax = subplots(figsize= isInside ? (2, 2) : (2, 2), layout="constrained")   #5,4
        styleAxis!(ax; xZoom=xLimitRecall, yRange=(0,yZoom), xStep=xStep, yStep=yStep, yScale=yScale)

        @info "Plotting PR curves"
        lastPlotData = nothing 
        for (idx, (legendLabel, currFilePR)) in enumerate(listFilePR)
            @info "Plot $idx: $legendLabel" file=currFilePR
            # If lineType is provided and nonempty, then pick its element if available.
            currentLineType = (length(lineTypes) ≥ idx && lineTypes[idx] != "") ? lineTypes[idx] : nothing
            currentLineWidth = (length(lineWidths) ≥ idx &&lineWidths[idx] != "") ? lineWidths[idx] : nothing
            lastPlotData = plotFileData!([ax], currFilePR, legendLabel, lineColors[idx]; 
                                 lineType=currentLineType, lineWidth = currentLineWidth, mode)
        end

        # Plot random PR line from the last file if available
        if lastPlotData !== nothing && lastPlotData[:randPR] !== nothing
            randPR = lastPlotData[:randPR]
            ax.axhline(randPR, linestyle="-.", linewidth=2, color=[0.6, 0.6, 0.6], label="Random")
        end

        # # Set labels and legend
        ax.set_xlabel("Recall", fontsize=axisTitleSize, family="Nimbus Sans")
        ax.set_ylabel("Precision", fontsize=axisTitleSize, family="Nimbus Sans")
        # save figure without legend
        if saveName !== nothing && !isempty(saveName)
            savePath = joinpath(dirOut, string(saveName, "_PR_noLegend.pdf"))
        else
            savePath = joinpath(dirOut, string(Dates.format(now(), "yyyymmdd_HHMMSS"), "_PR_noLegend.pdf"))
        end
        PyPlot.savefig(savePath, dpi=600)

        # Add Figure Legends
        if isInside
            ax.legend(borderaxespad=0.2, frameon=true)
        else
            ax.legend(loc="center", bbox_to_anchor=(1.25, 0.5), borderaxespad=0.2, frameon=false)
            fig.set_size_inches(3.5, 2.5)  # wider only for outside legend
        end
        savePathLegend = replace(savePath, "PR_noLegend.pdf" => "PR.pdf")
        PyPlot.savefig(savePathLegend, dpi=600)

    elseif length(yZoomPR) == 2
        # # ---- Two-subplot (broken y-axis) mode ----
        fig, (ax1, ax2) = subplots(2, 1, sharex=true, figsize=(2,2), gridspec_kw=Dict("height_ratios" => heightRatios,
                                "hspace" => 0.0), layout="constrained")   #6,5

        # Upper axis
        styleAxis!(ax1; xZoom=xLimitRecall, yRange=(yZoomPR[2],1.0), xStep=xStep, yStep=yStep, yScale=yScale, setXticks=false)
        # Lower axis
        styleAxis!(ax2; xZoom=xLimitRecall, yRange=(0,yZoomPR[1]), xStep=xStep, yStep=yStep, yScale=yScale)
        ax2.set_xlabel("Recall", fontsize=9)

        # Plot curves
        lastPlotData = nothing
        @info "Plotting PR curves"
        for (idx, (legendLabel, currFilePR)) in enumerate(listFilePR)
            @info "Plot $idx: $legendLabel; File: $currFilePR"
            # If lineType is provided and nonempty, then pick its element if available.
            currentLineType = (length(lineTypes) ≥ idx && lineTypes[idx] != "") ? lineTypes[idx] : nothing
            currentLineWidth = (length(lineWidths) ≥ idx &&lineWidths[idx] != "") ? lineWidths[idx] : nothing
            lastPlotData = plotFileData!([ax1, ax2], currFilePR, legendLabel, lineColors[idx]; 
                                 lineType=currentLineType, lineWidth = currentLineWidth, mode)

        end

        # Extract randPR from the last file
        # Plot the random PR line on both axes if available.
        if lastPlotData !== nothing && lastPlotData[:randPR] !== nothing
            randPR = lastPlotData[:randPR]
            # println(randPR)
            for ax in (ax1, ax2)
                ax.axhline(randPR, linestyle="-.", linewidth=2, color=[0.6, 0.6, 0.6], label="Random")
            end
        end

        # Hide the spines between ax1 and ax2
        ax1.spines["bottom"].set_visible(false)
        ax2.spines["top"].set_visible(false)
        ax1.tick_params(axis="x", which="both", bottom=false, top=false, labelbottom=false, labeltop=false)
        ax2.xaxis.tick_bottom()

        # Add break indicators
        d = 0.006  # Size of the diagonal lines in axes coordinates
        # Top-left and top-right diagonals for ax1
        ax1.plot([-d, +d], [-d, +d]; transform=ax1.transAxes, color="k", clip_on=false)  # Top-left diagonal
        ax1.plot([1 - d, 1 + d], [-d, +d]; transform=ax1.transAxes, color="k", clip_on=false)  # Top-right diagonal
        # Bottom-left and bottom-right diagonals for ax2
        ax2.plot([-d, +d], [1 - d, 1 + d]; transform=ax2.transAxes, color="k", clip_on=false)  # Bottom-left diagonal
        ax2.plot([1 - d, 1 + d], [1 - d, 1 + d]; transform=ax2.transAxes, color="k", clip_on=false)  # Bottom-right diagonal

        # Common labels: a shared y-axis label and x-axis label on the lower plot.
        # fig.text(0.01, 0.5, "Precision", va="center", rotation="vertical", fontsize=axisTitleSize)
        fig.supylabel("Precision", fontsize=axisTitleSize)
        ax2.set_xlabel("Recall", fontsize=axisTitleSize)


        # save figure without legend
        if saveName !== nothing && !isempty(saveName)
            savePath = joinpath(dirOut, string(saveName, "_PR_noLegend.pdf"))
        else
            savePath = joinpath(dirOut, string(Dates.format(now(), "yyyymmdd_HHMMSS"), "_PR_noLegend.pdf"))
        end
        PyPlot.savefig(savePath, dpi=600)

        # Add Figure Legends
        if isInside
            ax2.legend(fontsize= legendSize, borderaxespad=0.2, frameon=true)
        else
            ax2.legend(fontsize= legendSize, loc="center", bbox_to_anchor=(1.25, 0.5), borderaxespad=0.2, frameon=false)
            fig.set_size_inches(5, 4)  # wider only for outside legend
        end

        savePathLegend = replace(savePath, "PR_noLegend.pdf" => "PR.pdf")
        PyPlot.savefig(savePathLegend, dpi=600)

    end
    PyPlot.close("all")
end


# -------------------------------------------------------------------------------------
# Part 2: Making a DotPlot or BarPlot of AUPR using PyPlot 
# -------------------------------------------------------------------------------------
function plotAUPR(gsParam::OrderedDict{String,<:AbstractDict}, dirOut::String; 
                saveName::Union{Nothing, String} = nothing,
                metricType::String="full",
                figSize::Tuple{Real,Real}=(5,5),
                axisTitleSize::Int=9, tickLabelSize::Int=7,
                legendFontSize::Int=9, tickRotation::Int=45,
                plotType::String="",
                saveLegend::Bool=true)
    
    # ----------------------
    # Configure matplotlib
    # rc = PyPlot.matplotlib["rcParams"]
    # rc["font.family"] = "Nimbus Sans"
    # rc["axes.labelsize"] = axisTitleSize
    # rc["xtick.labelsize"] = tickLabelSize
    # rc["ytick.labelsize"] = tickLabelSize
    # rc["legend.fontsize"] = legendFontSize


    # Ensure output directory exists
    dirOut = joinpath(dirOut, "AUPR")
    mkpath(dirOut)

    # Load AUPR values
    function loadAupr(filePath::String, metricType::String)
        try
            data  = load(filePath)
            auprs = data["results"][:auprs]
            return lowercase(metricType) == "partial" ? auprs[:partial][:value] : auprs[:full]
        catch e
            @warn "Could not load AUPR from $filePath" exception=(e, catch_backtrace())
            return nothing
        end
    end
    
    lineColors = copy(DEFAULT_COLORS)
    # Build dataframe
    gsNames  = collect(keys(gsParam))
    netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
    
    dfRows = Vector{NamedTuple{(:xGroups,:Network,:AUPR),Tuple{String,String,Float64}}}()
    for gs in gsNames, net in netNames
        filePath = haskey(gsParam[gs], net) ? gsParam[gs][net] : nothing
        # And update the caller to skip nothing values
        auprVal = filePath === nothing ? nothing : loadAupr(filePath, metricType)
        push!(dfRows, (xGroups=gs, Network=net, AUPR=something(auprVal, NaN)))  # NaN is visible in plot
        # push!(dfRows, (xGroups=gs, Network=net, AUPR=auprVal))
    end
    
    df = DataFrame(dfRows)
    numGS  = length(unique(df.xGroups))
    numNet = length(unique(df.Network))
    
    # Decide plot type automatically if user hasn't overridden
    autoPlotType = ""
    if lowercase(plotType) in ["bar","dot"]
        autoPlotType = lowercase(plotType)
    else
        autoPlotType = (numGS == 1 || numNet == 1) ? "dot" : "bar"
    end
    
    xLabelVal = (numGS == 1 ? "Network" : "Gold Standards")
    yLabelVal = lowercase(metricType) == "partial" ? "Partial AUPR" : "AUPR"

    # Function to generate a plot
    function make_plot(includeLegend::Bool)
        fig, ax = plt.subplots(figsize=figSize, layout="constrained")
        colors = lineColors[1:numNet]
        
        if autoPlotType == "bar"
            if numGS == 1
                for (i, net) in enumerate(netNames)
                    idxs = findall(df.Network .== net)
                    ax.bar.(i, df.AUPR[idxs], color=colors[i], label=(includeLegend ? net : nothing))
                end
                ax.set_xticks(1:numNet)
                ax.set_xticklabels(netNames, rotation=tickRotation)
            elseif numNet == 1
                ax.bar.(1:numGS, df.AUPR, color=colors[1], label=(includeLegend ? netNames[1] : nothing))
                ax.set_xticks(1:numGS)
                ax.set_xticklabels(gsNames, rotation=tickRotation)
            else
                barWidth = 0.15
                groupWidth = numNet * barWidth
                groupSpacing = 0.4
                groupPositions = [i*(groupWidth + groupSpacing) for i in 0:(numGS-1)]
                for j in 1:numNet
                    offset = (j-(numNet+1)/2)*barWidth
                    xPos = [gp + groupWidth/2 + offset for gp in groupPositions]
                    auprs = [df.AUPR[(i-1)*numNet + j] for i in 1:numGS]
                    ax.bar(xPos, auprs, barWidth, color=colors[j], label=(includeLegend ? netNames[j] : nothing))
                end
                ax.set_xticks([gp + groupWidth/2 for gp in groupPositions])
                ax.set_xticklabels(gsNames, rotation=tickRotation)
            end
        elseif autoPlotType == "dot"
            if numGS == 1
                for (i, net) in enumerate(netNames)
                    idxs = findall(df.Network .== net)
                    ax.scatter.(i, df.AUPR[idxs], color=colors[i], s=40, label=(includeLegend ? net : nothing))
                end
                ax.set_xticks(1:numNet)
                ax.set_xticklabels(netNames, rotation=tickRotation)
            elseif numNet == 1
                ax.scatter.(1:numGS, df.AUPR, color=colors[1], s=40, label=(includeLegend ? netNames[1] : nothing))
                ax.set_xticks(1:numGS)
                
                ax.set_xticklabels(gsNames, rotation=tickRotation)
            else
                for j in 1:numNet
                    xvals = 1:numGS
                    yvals = [df.AUPR[(i-1)*numNet + j] for i in 1:numGS]
                    ax.scatter.(xvals, yvals, color=colors[j], s=40, label=(includeLegend ? netNames[j] : nothing))
                end
                ax.set_xticks(1:numGS)
                ax.set_xticklabels(gsNames, rotation=tickRotation)
            end
        end
        
        ax.set_xlabel(xLabelVal, fontsize=axisTitleSize)
        ax.set_ylabel(yLabelVal, fontsize=axisTitleSize)
        ax.tick_params(axis="both", which="both")
        # ax.tick_params(axis="both", which="both", labelsize=tickLabelSize, direction="in")
        ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
        if includeLegend
            # ax.legend(fontsize=legendFontSize, loc="best")
            ax.legend(fontsize=legendFontSize, loc="center", bbox_to_anchor=(1.25, 0.5), borderaxespad=0.2, frameon=false)
        else
            # Remove x-axis labels and ticks when no legend
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_xlabel("")
        end
        
        return fig, ax
    end
    
    saveName = (saveName === nothing || isempty(saveName)) ? "AUPR_$(Dates.format(now(), "yyyymmdd_HHMMSS"))" : saveName
    # Save plot according to saveLegend
    savePath = joinpath(dirOut, saveName * "_$(autoPlotType)_" * (saveLegend ? "withLegend" : "noLegend") * ".pdf")
    fig, _ = make_plot(saveLegend)
    fig.savefig(savePath, dpi=600)
    plt.close(fig)
    @info "Saved plot to $savePath"
end


# function plotAUPR(gsParam::OrderedDict{String,<:AbstractDict}, dirOut::String, saveName::String; figSize::Tuple{Real, Real})

#     """
#     Generate a visualization of Area Under the Precision-Recall Curve (AUPR) values using either a bar plot or a dot plot, depending on the input data structure.
#     - If there are multiple gold standards and networks, a bar plot is created, with each group of bars representing a gold standard and each bar within a group representing a network.
#     - If there is only one gold standard or one network, a dot plot is created, with each dot representing a network or gold standard.
#     - The plot is saved as a PDF file in the specified `dirOut` directory, with a resolution of 600 DPI.


#     # Arguments
#     - `gsParam::OrderedDict{String,Dict{String,String}}`: Mapping of gold standards to network file paths.
#     - `dirOut::String`: Directory to save the plot.
#     - `saveName::String`: Base name for the output file.

#     # Keywords
#     - `figSize::Tuple{Real, Real}=(5, 5)`: Size of the figure.

#     # Returns
#     Nothing, but saves the plot to a file.
#     """
        
#     function loadAupr(filePath)
#         try
#             data = load(filePath)
#             return data["results"][:auprs]
#         catch e
#             println("Error loading $filePath: $e")
#             return nothing
#         end
#     end

#     # Check if dirOut is either nothing or an empty string
#     dirOut = if dirOut === nothing || isempty(dirOut)
#         pwd()
#     else
#         mkpath(dirOut)
#     end

#     if isempty(figSize) || length(figSize) == 1
#         println("Warning: figSize must be a tuple of 2 element")
#         figSize = (5,5)
#     end

#     axisTitleSize = 16
#     tickLabelSize = 14
#     plotTitleSize = 16
#     legendSize = 10

#     lineColors = [ #Initial Color palette
#         "#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
#         "#984ea3", "#e41a1c", "#00ced1", "#000000",
#         "#5A9D5A", "#D96D3B", "#FFAD12", "#66628D", "#91569A",
#         "#B6742A", "#DD87B4", "#D26D7A", "#dede00"]

#     # Load aupr values
#     gsNames = collect(keys(gsParam))
#     netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
#     # auprValues = [loadAupr(gsParam[gs][netName]) !== nothing ? loadAupr(gsParam[gs][netName]) : 0.0 for gs in gsNames, netName in netNames]

#     auprValues = []
#     # Load AUPR values for each GS and file
#     for gs in gsNames
#         for netName in netNames
#             # Load the AUPR value
#             filePath = gsParam[gs][netName]
#             aupr = loadAupr(filePath)
            
#             # Use 0.0 if the AUPR value is missing
#             push!(auprValues, aupr === nothing ? 0.0 : aupr)
#         end
#     end
#     # Reshape auprValues into a matrix
#     auprMatrix = transpose(reshape(auprValues, length(netNames), length(gsNames)))

#     boolNet = any(length(netNames) > 1 for netNames in values(gsParam))
#     numNet = length(netNames)
#     numGS = length(gsNames)
#     # =
#     # ---- Load AUPR values (Works but not using)
#     # auprData = Dict(gs => Dict(name => loadAupr(path) for (name, path) in files) for (gs, files) in gsParam)
#     # netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
#     # auprValues = [get(auprData[gs], netName, 0.0) for gs in gsNames, netName in netNames]
#     = #
#     if boolNet && numGS > 1
#         # Create a figure
#         fig, ax = plt.subplots(figsize=figSize, layout="constrained")
#         # Define parameters for grouping:
#         barWidth = 0.15               # Width of each bar.
#         groupSpacing = 0.4            # Extra space between groups.
#         groupWidth = numNet * barWidth # Total width occupied by bars in one group.

#         # Compute positions for each group.
#         # Compute a vector where each element is the left edge (or starting point) of each group.
#         # Use 0-based indexing for groups.
#         groupPositions = [i * (groupWidth + groupSpacing) for i in 0:(numGS-1)]

#         # Now plot each bar (each network) in every group.
#         # Use a centered offset for each bar in the group.
#         for idx in 1:numNet
#             # Calculate the offset to center the bar within the group.
#             # (idx - (numNet+1)/2)*barWidth shifts bars so they cluster around the center.
#             offset = (idx - (numNet+1)/2)*barWidth
#             # The x-positions for the bars are:
#             #   groupPositions + offset + groupWidth/2
#             # Adding groupWidth/2 centers the bars in each group.
#             x_positions = [gp + groupWidth/2 + offset for gp in groupPositions]
#             # Plot the bar for curr network across all groups.
#             ax.bar(x_positions, auprMatrix[:, idx], barWidth, label=netNames[idx],
#                 color = lineColors[mod1(idx, length(lineColors))])
#         end

#         # Set the xticks to the center positions of each group.
#         groupCenters = [gp + groupWidth/2 for gp in groupPositions]
#         ax.set_xticks(groupCenters)
#         ax.set_xticklabels(gsNames)
#         ax.set_xlabel("Gold Standards", fontsize=axisTitleSize)
#         ax.set_ylabel("AUPR", fontsize=axisTitleSize)
#         ax.tick_params(axis="both", which="both", labelsize=tickLabelSize)
#         ax.legend(fontsize=legendSize, loc="center left", bbox_to_anchor=(1, 0.5))

#     elseif numGS == 1 || numNet == 1
#         # Dot Plot
#         fig, ax = subplots(figsize= figSize, layout="constrained")
#         if numGS == 1
#             # One GS, multiple files
#             for idx in 1:numNet
#                 ax.scatter(idx, auprMatrix[idx], color=lineColors[idx], label=netNames[idx])
#                 ax.text(idx, auprMatrix[idx], netNames[idx], fontsize=9, ha="right")
#             end
#             ax.set_xlabel("Network", fontsize=axisTitleSize)
#             ax.set_xticks([])
#             ax.set_xticklabels([])
#             # ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#         else
#             for idx in 1:numGS
#                 ax.scatter(idx, auprMatrix[idx], color=lineColors[idx], label=gsNames[idx])
#                 ax.text(idx, auprMatrix[idx], gsNames[idx], fontsize=9, ha="right")
#             end
#             ax.set_xlabel("Gold Standard", fontsize=axisTitleSize)
#             ax.set_xticks([])
#             ax.set_xticklabels([])
#             # ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#         end
#     end 

#     # Common settings for both plots.
#     ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#     ax.tick_params(axis="both", which="both", labelsize=tickLabelSize)
#     # plt.tight_layout()


#     if !isempty(saveName)
#         plt.savefig(joinpath(dirOut, saveName * "_AUPR.pdf"), dpi=600)
#     else
#         dateStr = Dates.format(now(), "yyyymmdd_HHMMSS")
#         plt..savefig(joinpath(dirOut, dateStr * "_AUPR.pdf"), dpi=600)
#     end

#     plt.close("all")
# end






# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------ Using StatsPlots
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# function plotAUPR(gsParam::OrderedDict{String,<:AbstractDict}, dirOut::String, 
#     saveName::String, plotType::String = "bar"; 
#     figSize::Tuple{Real,Real} = (8,5), axisTitleSize::Int = 13,
#     tickLabelSize::Int = 10)

#     legendFontSize = 9

#     # Check if dirOut is either nothing or an empty string
#     dirOut = if dirOut === nothing || isempty(dirOut)
#         pwd()
#     else
#         mkpath(dirOut)
#     end

#     # Dummy load function – replace with your actual data-loading mechanism.
#     function loadAupr(filePath)
#         try
#             data = load(filePath)
#             return data["results"][:auprs]
#         catch e
#             println("Error loading $filePath: $e")
#             return nothing
#         end
#     end

#     # Initial color palette
#     lineColors = ["#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
#             "#984ea3", "#e41a1c", "#00ced1", "#000000",
#             "#5A9D5A", "#D96D3B", "#FFAD12", "#66628D", "#91569A",
#             "#B6742A", "#DD87B4", "#D26D7A", "#dede00"]

#     # Build a DataFrame with columns: xGroups, Network, AUPR
#     gsNames = collect(keys(gsParam))
#     netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
#     if length(netNames) > length(lineColors)
#         lineColors = padColors(netNames)  # Ensure padColors is defined!
#     end

#     dfRows = Vector{NamedTuple{(:xGroups, :Network, :AUPR), Tuple{String,String,Float64}}}()
#     for gs in gsNames
#         for net in netNames
#             filePath = haskey(gsParam[gs], net) ? gsParam[gs][net] : nothing
#             temp = filePath === nothing ? nothing : loadAupr(filePath)
#             val = (temp === nothing || temp === missing) ? 0.0 : temp
#             push!(dfRows, (xGroups = gs, Network = net, AUPR = val))
#         end
#     end
#     df = DataFrame(dfRows)

#     # Number of unique groups.
#     numGS  = length(unique(df.xGroups))
#     numNet = length(unique(df.Network))

#     # Set common labels and legend; use camelCase.
#     xLabelVal = (numGS == 1 ? "Network" : "Gold Standards")
#     yLabelVal = "AUPR"
#     legendPos = :outerright

#     # Common arguments (for scatter and bar, without color arguments).
#     commonArgs = (xlabel = xLabelVal, ylabel = yLabelVal, legend = legendPos, 
#             guidefontsize = axisTitleSize, tickfontsize = tickLabelSize, legendfontsize = legendFontSize,
#             framestyle = :box)

#     colors = lineColors[1:numNet]

#     p = nothing
#     if lowercase(plotType) == "scatter"
#         if numGS == 1
#             p = @df df sp.scatter(:Network, :AUPR; markersize = 6, 
#                 color = colors, xrotation = 45, commonArgs...)
#         elseif numNet == 1
#             p = @df df sp.scatter(:xGroups, :AUPR; markersize = 6,
#                 color = lineColors[1], xrotation = 45, commonArgs...)
#         else
#             p = @df df sp.scatter(:xGroups, :AUPR; group = :Network, markersize = 6,
#                 palette = colors, commonArgs...)
#         end

#     elseif lowercase(plotType) == "bar"
#         if numGS == 1
#             p = @df df sp.bar(:Network, :AUPR; 
#                 color = colors, xrotation = 45, commonArgs...)
#         elseif numNet == 1
#             p = @df df sp.bar(:xGroups, :AUPR; 
#                 color = lineColors[1], xrotation = 45, commonArgs...)
#         else
#             p = @df df sp.groupedbar(:xGroups, :AUPR; group = :Network,
#                 palette = colors, commonArgs...)
#         end
#     else
#         error("Invalid plot type. Please choose 'scatter' or 'bar'")
#     end

#     # Adjust figure size (e.g. figSize = (6,3.5) converts to (600,350) pixels)
#     p = sp.plot(p, size = (Int(figSize[1]*100), Int(figSize[2]*100)), titlefontsize = axisTitleSize)


#     if !isempty(saveName)
#         savePath = joinpath(dirOut, saveName * "_AUPR_$(plotType).pdf")
#         sp.savefig(p, savePath)
#     else
#         dateStr = Dates.format(now(), "yyyymmdd_HHMMSS")
#         savePath = joinpath(dirOut, dateStr * "_AUPR_$(plotType).pdf")
#         sp.savefig(p, savePath)
#     end
# end



