# =============================================================================
#  dataUtils.jl — Examples for InferelatorJL data utility functions
#
#  What:
#    Demonstrates data manipulation utilities using small synthetic datasets.
#    No real input files required — all data is generated inline.
#    Sections:
#      1. DataUtils         — reshape, normalize, binarize, merge prior matrices
#      2. PartialCorrelation — precision-matrix and regression-based methods
#
#  Required inputs:  None (all data is synthetic)
#
#  Expected outputs: Printed results to the REPL / stdout.
#    Section 1 also writes a temporary file to tempdir().
#
#  Usage:
#    julia examples/dataUtils.jl
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

# ── matrixToEdgeList / edgeListToMatrix ──────────────────────────────────────
# Wide prior matrix: genes (rows) × TFs (columns)
widePrior = DataFrame(
    Gene = ["Gata3", "Foxp3", "Tbx21", "Rorc"],
    TF_A = [1.0, 0.0, 0.5, 0.0],
    TF_B = [0.0, 1.0, 0.0, 1.0],
    TF_C = [0.8, 0.0, 0.0, 0.6]
)

# Convert wide matrix to edge list (TF, Gene, Weight)
longPrior = matrixToEdgeList(widePrior)
println("Wide matrix → Edge list (first 5 rows):")
println(first(longPrior, 5))

# Round-trip back to wide matrix
widePrior2 = edgeListToMatrix(longPrior; indices = (2, 1, 3))
println("\nEdge list → Wide matrix (restored):")
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
sparseFile = joinpath(tmpDir, "prior_example_sp.tsv")   # _sp = long/edge list format
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
