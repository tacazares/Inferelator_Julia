module Metrics

using JLD2, InlineStrings
using StatsPlots
using PyPlot
using Colors
using Dates
using DataFrames
using OrderedCollections
using Measures
using CSV
using Base.Threads
using Interpolations
using Statistics

# ----------------------------
# Plot defaults
# ----------------------------

const sp = StatsPlots
const dpi = 600

function setPlotDefaults!()
    sp.default(
        fontfamily = "Nimbus Sans",
        guidefont  = font(9, "Nimbus Sans"),
        tickfont   = font(7, "Nimbus Sans"),
        legendfont = font(9, "Nimbus Sans"),
        dpi        = dpi,
        margin     = 5
    )

    rc = PyPlot.matplotlib.rcParams

    rc["font.family"]      = "Nimbus Sans"
    rc["axes.titlesize"]   = 9
    rc["axes.labelsize"]  = 9
    rc["xtick.labelsize"] = 7
    rc["ytick.labelsize"] = 7
    rc["legend.fontsize"] = 9

    rc["figure.dpi"]  = dpi
    rc["savefig.dpi"] = dpi

    return nothing
end


# ----------------------------
# Public API
# ----------------------------

export setPlotDefaults!,
       plotPRCurve,
       plotROCCurve,
       plotPRCurves,
       plotAUPR,
       computeMacroMetrics,
       evaluatePerTF,
       computePR


# ----------------------------
# Sub-files
# ----------------------------
include("constants.jl")  
include("utils.jl")  
include("calcPRinfTRNs.jl")
include("plotting/plotSingleUtils.jl")
include("plotting/plotBatchMetrics.jl")


end # module Metrics