using JLD2, InlineStrings
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

const dpi = 600

function setPlotDefaults!()
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