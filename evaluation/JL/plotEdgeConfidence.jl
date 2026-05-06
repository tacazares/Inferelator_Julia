using CSV, DataFrames, OrderedCollections, StatsBase, PyPlot, Dates

# Set global font defaults
PyPlot.matplotlib.rcParams["font.family"]       = "Nimbus Sans"
PyPlot.matplotlib.rcParams["axes.titlesize"]    = 12
PyPlot.matplotlib.rcParams["axes.labelsize"]    = 10
PyPlot.matplotlib.rcParams["xtick.labelsize"]   = 8
PyPlot.matplotlib.rcParams["ytick.labelsize"]   = 8
PyPlot.matplotlib.rcParams["legend.fontsize"]   = 9

# ----------------------------------------------------------------------
# Helper: load and process stability values from an edges file
# Returns (normVal, nbinsUse) or (nothing, nothing) if file missing/invalid
#
# Expects bStARS or bEBIC (count-scale) edges files: Stability = integer
# subsample count + fractional correlation weight. floor() recovers the count;
# one bin per integer value (bin width = 1 subsample).
# nbins overrides the automatic bin count when supplied.
# ----------------------------------------------------------------------
function _loadConfidence(file::String; normalize::Bool=false, nbins::Union{Nothing,Int}=nothing)
    if !isfile(file)
        @warn "File not found: $file"
        return nothing, nothing
    end
    df = try
        CSV.read(file, DataFrame; delim='\t')
    catch e
        @warn "Error reading $file" exception=e
        return nothing, nothing
    end
    if !hasproperty(df, :Stability)
        @warn "No 'Stability' column in $file"
        return nothing, nothing
    end
    vals    = Float64.(df[:, :Stability])
    floored = floor.(vals)
    maxVal  = maximum(floored)
    minVal  = minimum(floored)

    normVal  = normalize ? (maxVal == 0 ? floored : floored ./ maxVal) : floored
    nbinsUse = isnothing(nbins) ? Int(minVal == 0 ? maxVal + 1 : maxVal) : nbins

    return normVal, nbinsUse
end

# ----------------------------------------------------------------------
# histogramConfidencesStacked
#
# Creates confidence histograms from a named dict of edge files.
#
# Modes:
#   layered=false  →  one PDF per network, saved to dirOut
#   layered=true   →  one overlaid PDF for all networks, saved to dirOut
#
# Arguments:
#   netFiles   : OrderedDict mapping label => path to edges.tsv
#   dirOut     : output directory (created if needed)
#   saveName   : base filename (timestamp used if nothing)
#   nbins      : override bin count (auto if nothing)
#   layered    : true = overlaid plot, false = individual plots
#   normalize  : normalize stability values to [0,1]
#   logscale   : log10 y-axis
# ----------------------------------------------------------------------
function histogramConfidencesStacked(netFiles::OrderedDict{String,String},
                                     dirOut::String;
                                     saveName::Union{Nothing,String}=nothing,
                                     nbins::Union{Nothing,Int}=nothing,
                                     layered::Bool=true,
                                     normalize::Bool=false,
                                     logscale::Bool=false)

    dirOut = (isnothing(dirOut) || isempty(dirOut)) ? pwd() : mkpath(dirOut)

    if !layered
        # ---- Individual plots: one PDF per network ----
        for (netName, file) in netFiles
            normVal, nbinsUse = _loadConfidence(file; normalize, nbins)
            isnothing(normVal) && continue

            fig, ax = subplots(figsize=(5, 4), layout="constrained")
            ax.hist(normVal, bins=nbinsUse, color="steelblue", edgecolor="black", alpha=0.9)
            ax.set_title(netName)
            ax.set_xlabel(normalize ? "Normalized Confidence" : "Confidence (Subsample Count)")
            ax.set_ylabel("Number of TF-Gene edges")
            if logscale
                ax.set_yscale("log")
            end
            ax.minorticks_on()
            ax.grid(true, which="major", linestyle="-",  linewidth=0.5, color="lightgray")
            ax.grid(true, which="minor", linestyle=":",  linewidth=0.2, color="#e8e8e8")

            suffix   = logscale ? "_logscale" : ""
            outPath  = joinpath(dirOut, netName * "_hist_confidence_" * string(nbinsUse) * suffix * ".pdf")
            PyPlot.savefig(outPath, dpi=600)
            PyPlot.close(fig)
            println("Saved: $outPath")
        end

    else
        # ---- Layered plot: all networks overlaid ----
        allVals      = Dict{String, Vector{Float64}}()
        globalMaxRaw = 0.0

        for (netName, file) in netFiles
            normVal, _ = _loadConfidence(file; normalize, nbins)
            isnothing(normVal) && continue
            allVals[netName] = normVal
            globalMaxRaw = max(globalMaxRaw, maximum(floor.(normVal)))
        end

        if isempty(allVals)
            @warn "No data found in any file."
            return
        end

        nbinsUse = isnothing(nbins) ? (globalMaxRaw > 1 ? Int(globalMaxRaw) + 1 : 30) : nbins
        println("Using $nbinsUse bins (max raw confidence: $globalMaxRaw).")

        fig, ax = subplots(figsize=(5, 4), layout="constrained")
        for (netName, vals) in allVals
            ax.hist(vals, bins=nbinsUse, alpha=0.6, label=netName, edgecolor="black", linewidth=0.3)
        end
        ax.set_xlabel(normalize ? "Normalized Confidence" : "Confidence (Subsample Count)")
        ax.set_ylabel("Number of TF-Gene edges")
        if logscale
            ax.set_yscale("log")
        end
        ax.legend(frameon=false)
        ax.minorticks_on()
        ax.grid(true, which="major", linestyle="-",  linewidth=0.5, color="lightgray")
        ax.grid(true, which="minor", linestyle=":",  linewidth=0.3, color="lightgray")

        ts      = Dates.format(now(), "yyyymmdd_HHMMSS")
        name    = (isnothing(saveName) || isempty(saveName)) ? ts : saveName
        suffix  = logscale ? "_logscale" : ""
        outPath = joinpath(dirOut, name * "_Confidences" * suffix * ".pdf")
        PyPlot.savefig(outPath, dpi=600)
        PyPlot.close(fig)
        println("Saved: $outPath")
    end
end

# ----------------------------------------------------------------------
# histogramConfidencesDir
#
# Generates individual confidence histograms for TFA and TFmRNA
# subfolders found inside each directory in currNetDirs.
# Plots are saved inside each subfolder.
#
# Arguments:
#   currNetDirs : vector of network output directories
#   normalize   : normalize stability values to [0,1]
#   logscale    : log10 y-axis
# ----------------------------------------------------------------------
function histogramConfidencesDir(currNetDirs::Vector{String};
                                 normalize::Bool=false,
                                 logscale::Bool=false,
                                 nbins::Union{Nothing,Int}=nothing)

    for currNetDir in currNetDirs
        subfolders    = filter(isdir, readdir(currNetDir; join=true))
        targetFolders = filter(s -> basename(s) in ["TFA", "TFmRNA"], subfolders)

        for subfolder in targetFolders
            filePath = joinpath(subfolder, "edges.tsv")
            normVal, nbinsUse = _loadConfidence(filePath; normalize, nbins)
            isnothing(normVal) && continue

            fig, ax = subplots(figsize=(5, 4), layout="constrained")
            ax.hist(normVal, bins=nbinsUse, color="steelblue", edgecolor="black", alpha=0.9)
            ax.set_title(basename(subfolder))
            ax.set_xlabel(normalize ? "Normalized Confidence" : "Confidence (Subsample Count)")
            ax.set_ylabel("Number of TF-Gene edges")
            if logscale
                ax.set_yscale("log")
            end
            ax.minorticks_on()
            ax.grid(true, which="major", linestyle="-",  linewidth=0.5, color="lightgray")
            ax.grid(true, which="minor", linestyle=":",  linewidth=0.2, color="#e8e8e8")

            suffix   = logscale ? "_logscale" : ""
            saveFile = joinpath(subfolder, "confidence_distribution_" * string(nbinsUse) * suffix * ".pdf")
            PyPlot.savefig(saveFile, dpi=600)
            PyPlot.close(fig)
            println("Saved: $saveFile")
        end
    end
end


# ======================================================================
# USAGE
# ======================================================================

# --- histogramConfidencesDir: one plot per TFA/TFmRNA subfolder ---
# currNetDirs = [
#     "/path/to/network/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63",
#     "/path/to/network/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63"
# ]
# histogramConfidencesDir(currNetDirs; logscale=false, normalize=false)


# --- histogramConfidencesStacked: individual plots per network ---
# netFiles = OrderedDict(
#     "Network A" => "/path/to/networkA/TFA/edges.tsv",
#     "Network B" => "/path/to/networkB/TFA/edges.tsv"
# )
# histogramConfidencesStacked(netFiles, "/path/to/output";
#                              saveName="myNetworks", layered=false,
#                              normalize=false, logscale=false)


# --- histogramConfidencesStacked: overlaid plot for all networks ---
# histogramConfidencesStacked(netFiles, "/path/to/output";
#                              saveName="myNetworks_overlaid", layered=true,
#                              normalize=false, logscale=false)
