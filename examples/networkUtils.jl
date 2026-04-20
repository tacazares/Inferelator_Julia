# =============================================================================
#  networkUtils.jl — Examples for InferelatorJL network utility functions
#
#  What:
#    Demonstrates network I/O and aggregation utilities using synthetic edge data.
#    No real input files required — all data is generated inline.
#    Sections:
#      1. NetworkIO         — write edge tables and save data structs
#      2. aggregateNetworks — combine multiple GRN edge files into one network
#
#  Required inputs:  None (all data is synthetic)
#
#  Expected outputs: Printed results to the REPL / stdout.
#    Both sections write temporary files to tempdir().
#
#  Usage:
#    julia examples/networkUtils.jl
#    or paste individual sections into a Julia REPL
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================

using InferelatorJL
using DataFrames, CSV


# =============================================================================
# 1. NetworkIO
# =============================================================================

println("\n===== 1. NetworkIO =====\n")

# ── writeNetworkTable! ────────────────────────────────────────────────────────
# Populate a minimal BuildGrn with synthetic edge data and write TSV outputs.
buildGrn = BuildGrn()
buildGrn.regs               = ["TF_A", "TF_B", "TF_A", "TF_C"]
buildGrn.targs              = ["Gata3", "Foxp3", "Tbx21", "Rorc"]
buildGrn.signedQuantile     = [0.92, -0.75, 0.60, 0.85]
buildGrn.rankings           = [0.88, 0.72, 0.55, 0.80]
buildGrn.partialCorrelation = [0.45, -0.38, 0.30, 0.41]
buildGrn.inPrior            = ["Yes", "Yes", "No", "Yes"]
buildGrn.networkMat         = hcat(buildGrn.regs, buildGrn.targs,
                                    buildGrn.signedQuantile, buildGrn.rankings,
                                    buildGrn.partialCorrelation, buildGrn.inPrior)
buildGrn.networkMatSubset   = buildGrn.networkMat[1:3, :]  # top 3 edges as subset

tmpDir = tempdir()
netDir = joinpath(tmpDir, "network_example")
mkpath(netDir)
writeNetworkTable!(buildGrn; outputDir = netDir)
println("edges.tsv written to: $netDir")
println(read(joinpath(netDir, "edges.tsv"), String))


# ── saveData ──────────────────────────────────────────────────────────────────
# Save all four core structs to a .jld2 file for checkpointing.
# Reload in a fresh session with:
#   using InferelatorJL, JLD2
#   @load joinpath(netDir, "checkpoint.jld2") expressionData tfaData grnData buildGrn
#
# (Skipped here since expressionData/tfaData/grnData require real inputs)


# =============================================================================
# 2. aggregateNetworks
# =============================================================================

println("\n===== 2. aggregateNetworks =====\n")

# Write two synthetic edge files (TFA and mRNA modes) then combine them.
tfaEdges = DataFrame(
    TF             = ["TF_A", "TF_B", "TF_C", "TF_A"],
    Gene           = ["Gata3", "Foxp3", "Tbx21", "Rorc"],
    signedQuantile = [0.92, -0.75, 0.60, 0.85],
    Stability      = [0.88, 0.72, 0.55, 0.80],
    Correlation    = [0.45, -0.38, 0.30, 0.41],
    inPrior        = ["Yes", "Yes", "No", "Yes"]
)

mrnaEdges = DataFrame(
    TF             = ["TF_A", "TF_B", "TF_C", "TF_D"],
    Gene           = ["Gata3", "Foxp3", "Tbx21", "Rorc"],
    signedQuantile = [0.80, -0.65, 0.70, 0.55],
    Stability      = [0.75, 0.60, 0.68, 0.50],
    Correlation    = [0.40, -0.32, 0.35, 0.28],
    inPrior        = ["Yes", "Yes", "No", "No"]
)

tfaFile  = joinpath(tmpDir, "TFA_edges.tsv")
mrnaFile = joinpath(tmpDir, "mRNA_edges.tsv")
CSV.write(tfaFile,  tfaEdges;  delim = '\t')
CSV.write(mrnaFile, mrnaEdges; delim = '\t')

combDir = joinpath(tmpDir, "Combined")

# Combine using max stability per (TF, Gene) pair
combinedMax = aggregateNetworks(
    [tfaFile, mrnaFile];
    method              = :max,
    meanEdgesPerGene    = 3,
    useMeanEdgesPerGene = true,
    outputDir           = combDir
)
println("Combined network (:max strategy), $(nrow(combinedMax)) edges:")
println(combinedMax)

# Combine using mean stability
combinedMean = aggregateNetworks(
    [tfaFile, mrnaFile];
    method              = :mean,
    meanEdgesPerGene    = 3,
    useMeanEdgesPerGene = true,
    outputDir           = combDir
)
println("\nCombined network (:mean strategy), $(nrow(combinedMean)) edges:")
println(combinedMean)
