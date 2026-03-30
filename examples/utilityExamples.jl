# =============================================================================
#  utilityExamples.jl — Self-contained examples for InferelatorJL utility functions
#
#  What:
#    Demonstrates all exported utility functions using small synthetic datasets.
#    No real input files are required — all data is generated inline.
#    Sections:
#      1. DataUtils         — reshape, normalize, binarize, merge prior matrices
#      2. PartialCorrelation — precision-matrix and regression-based methods
#      3. NetworkIO         — write edge tables and save data structs
#      4. aggregateNetworks  — combine multiple GRN edge files into one network
#
#  Required inputs:  None (all data is synthetic)
#
#  Expected outputs: Printed results to the REPL / stdout.
#    Sections 3 and 4 also write temporary files to tempdir().
#
#  Usage:
#    julia examples/utilityExamples.jl
#    or paste individual sections into a Julia REPL
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================

using InferelatorJL
using DataFrames, CSV, Random, Statistics, LinearAlgebra


# =============================================================================
# 1. DataUtils
# =============================================================================

println("\n===== 1. DataUtils =====\n")

# ── convertToLong / convertToWide ────────────────────────────────────────────
# Wide prior matrix: genes (rows) × TFs (columns)
widePrior = DataFrame(
    Gene = ["Gata3", "Foxp3", "Tbx21", "Rorc"],
    TF_A = [1.0, 0.0, 0.5, 0.0],
    TF_B = [0.0, 1.0, 0.0, 1.0],
    TF_C = [0.8, 0.0, 0.0, 0.6]
)

# Convert to long (TF, Gene, Weight)
longPrior = convertToLong(widePrior)
println("Wide → Long (first 5 rows):")
println(first(longPrior, 5))

# Round-trip back to wide
widePrior2 = convertToWide(longPrior; indices = (2, 1, 3))
println("\nLong → Wide (restored):")
println(widePrior2)


# ── frobeniusNormalize ────────────────────────────────────────────────────────
# Normalize each column of the prior matrix so each column has unit L2 norm.
# Typical use: normalize prior before feeding into penalty matrix construction.
priorNorm = frobeniusNormalize(widePrior, :column)
println("\nColumn-normalized prior (each TF column has ||·||₂ = 1):")
println(priorNorm)

# Verify norms are 1
details, nTrue, nFalse = check_column_norms(priorNorm; atol = 1e-8)
println("Columns with unit norm: $nTrue / $(nTrue + nFalse)")

# Row normalization (normalize each gene's regulatory weight vector)
priorNormRow = frobeniusNormalize(widePrior, :row)
println("\nRow-normalized prior:")
println(priorNormRow)


# ── binarizeNumeric! ──────────────────────────────────────────────────────────
# Replace all non-zero values with 1 (convert continuous prior → binary prior).
priorBinary = deepcopy(widePrior)
binarizeNumeric!(priorBinary)
println("\nBinarized prior:")
println(priorBinary)


# ── mergeDFs ──────────────────────────────────────────────────────────────────
# Merge two prior DataFrames from different assays (ATAC + ChIP) by summing.
# Useful when combining evidence from multiple prior sources.
prior_atac = DataFrame(Gene = ["Gata3","Foxp3","Tbx21"], TF_A = [1.0,0.0,0.5], TF_B = [0.0,1.0,0.0])
prior_chip = DataFrame(Gene = ["Foxp3","Tbx21","Rorc"],  TF_A = [0.0,0.3,0.0], TF_C = [1.0,0.0,1.0])

mergedPrior = mergeDFs([prior_atac, prior_chip], :Gene, "sum")
println("\nMerged prior (ATAC + ChIP, sum):")
println(mergedPrior)

mergedPriorAvg = mergeDFs([prior_atac, prior_chip], :Gene, "avg")
println("\nMerged prior (ATAC + ChIP, average):")
println(mergedPriorAvg)


# ── completeDF ────────────────────────────────────────────────────────────────
# Align a DataFrame to a fixed set of row IDs and column names (fills missing → 0).
allGenes = ["Gata3", "Foxp3", "Tbx21", "Rorc", "Il2ra"]
allTFs   = [:TF_A, :TF_B, :TF_C, :TF_D]

completedPrior = completeDF(prior_atac, :Gene, allGenes, allTFs)
println("\ncompleteDF — prior_atac aligned to full gene/TF universe:")
println(completedPrior)


# ── writeTSVWithEmptyFirstHeader ──────────────────────────────────────────────
# Write a DataFrame to TSV with the first (row-label) column header left blank.
# This is the sparse prior format expected by processPriorFile!.
tmpDir     = tempdir()
sparseFile = joinpath(tmpDir, "prior_example_sp.tsv")
writeTSVWithEmptyFirstHeader(priorBinary, sparseFile; delim = '\t')
println("\nSparse prior written to: $sparseFile")
println(read(sparseFile, String)[1:min(200, filesize(sparseFile))])


# =============================================================================
# 2. PartialCorrelation
# =============================================================================

println("\n===== 2. PartialCorrelation =====\n")

# Synthetic expression matrix: 50 cells × 6 genes
Random.seed!(42)
X = randn(50, 6)
# Add a latent signal so genes 1 and 3 are partially correlated
latent = randn(50)
X[:, 1] .+= 0.8 .* latent
X[:, 3] .+= 0.7 .* latent

# ── Method 1: Precision matrix ────────────────────────────────────────────────
# partialCorrelationMat returns the full p×p partial correlation matrix.
# Note: accessed via module prefix since it is an internal function.
Pfull = InferelatorJL.partialCorrelationMat(X; epsilon = 1e-6, first_vs_all = false)
println("Full partial correlation matrix (6×6):")
display(round.(Pfull, digits = 3))

# first_vs_all = true: only partial correlations of column 1 vs all others.
# This is how it is used inside rankEdges! (target gene vs. TF predictors).
P1 = InferelatorJL.partialCorrelationMat(X; epsilon = 1e-6, first_vs_all = true)
println("\nPartial correlations of gene 1 vs. all others (1×6):")
println(round.(P1, digits = 3))

# ── Method 2: Regression residuals ───────────────────────────────────────────
# partialCorrReg computes the same quantity via OLS residuals.
# Slower than the precision-matrix method; use for small matrices.
Preg = InferelatorJL.partialCorrReg(X; first_vs_all = false)
println("\nPartial correlation (regression method, 6×6):")
display(round.(Preg, digits = 3))

println("\nDifference between methods (should be ~0):")
println("Max abs diff: ", round(maximum(abs.(Pfull[2:end, 2:end] .- Preg[2:end, 2:end])), digits = 6))


# =============================================================================
# 3. NetworkIO
# =============================================================================

println("\n===== 3. NetworkIO =====\n")

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
# (Skipped here since data/tfaData/grnData require real inputs)


# =============================================================================
# 4. aggregateNetworks
# =============================================================================

println("\n===== 4. aggregateNetworks =====\n")

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
