# src/data/PseudobulkData.jl
#
# Pseudobulk aggregation, normalization, and batch correction utilities.
# Used to prepare single-cell expression data as input to InferelatorJL.
#
# `pseudoBulk`  accepts any AnnData-like object (e.g. from Muon.jl) via
# duck typing — Muon.jl is not a package dependency but must be loaded by
# the user before calling this function.
#
# RCall-based functions (vstNormalizeR, removeBatchEffectR) are NOT included
# here to avoid a heavy R dependency. Use the workflow scripts in
# examples/prepareData/pseudobulkWorkflow.jl for the exact DESeq2 + limma pipeline.


# ============================================================
#  PSEUDOBULK AGGREGATION
# ============================================================

"""
    pseudoBulk(adata, identVar::String) -> DataFrame

Aggregate single-cell counts into pseudobulk by summing across all cells
in each group defined by `identVar`.

`adata` is any AnnData-like object with:
- `.X`   : cells × features sparse count matrix
- `.obs` : cell metadata DataFrame (must contain `identVar` column)
- `.var` : feature metadata DataFrame (must contain a `"name"` column)

Typically loaded via `Muon.readh5ad`. All-zero rows are removed.

# Arguments
- `adata`    : AnnData object loaded via Muon.jl
- `identVar` : column in `adata.obs` to group cells by (e.g. `"celltype_rep"`)

# Returns
DataFrame with a leading `:RowNames` column and one column per unique group.

# Example
```julia
using Muon, InferelatorJL
adata  = readh5ad("data.h5ad")
pbRaw  = pseudoBulk(adata, "celltype_rep")
```
"""


# reference: https://docs.juliahub.com/Muon/QfqCh/0.2.1/

function createPseudobulk(adata, identVar::String, saveDir::String)
    # Read the AnnData object
    #adata = h5open(dataset, "r")
    # adata = read(adata)

    spMat = adata.X  #sparse matrix (cells by feature)
    
    identCol = adata.obs[!, identVar]  # element of column corresponding to identVar. Note: adata.obs -> metadata::dataframe
    uniqueIdents = unique(identCol)

    features = adata.var[!, "name"]  # contains rownames
    pseudobulkDF = DataFrame(zeros(length(features), length(uniqueIdents)), uniqueIdents)

    for ident in uniqueIdents
        currIndices = findall(x -> x == ident, identCol)
        # Subset data for the current identity using logical indexing
        currData = spMat[currIndices, :]
        
        # sum counts across all cells per row/feature
        currCounts = sum(currData, dims=1)
        pseudobulkDF[!, ident] .= vec(currCounts)
    end

    mat = Matrix(pseudobulkDF);
    sumsRow = vec(sum(mat, dims=2));  # Sum along rows;

    pseudobulkDF = insertcols!(pseudobulkDF, 1, :RowNames => features) 
    pseudobulkDF = pseudobulkDF[sumsRow .> 0, :]

    # Save the pseudobulk DataFrame
    currTime = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    outPath = joinpath(saveDir, "PseudobulkData_$(identVar)_$(currTime).tsv")
    # open(outPath, "w") do io
    #     write(io, pseudobulkDF)
    # end
    CSV.write(outPath, string.(pseudobulkDF); delim = '\t')

    println("Pseudobulk DataFrame saved to: ", outPath)
    return pseudobulkDF 
end


function pseudoBulk(adata, identVar::String)
    spMat        = adata.X
    identCol     = adata.obs[!, identVar]
    uniqueIdents = unique(identCol)
    features     = adata.var[!, "name"]

    pb = DataFrame(:RowNames => features)
    for ident in uniqueIdents
        idx    = findall(==(ident), identCol)
        currData = spMat[idx, :]
        counts = vec(sum(currData, dims=1))
        # counts = vec(sum(spMat[idx, :], dims=1))
        pb[!, Symbol(ident)] = Float64.(counts)
    end

    mat     = Matrix(pb[:, 2:end])
    rowSums = vec(sum(mat, dims=2))
    return pb[rowSums .> 0, :]
end


# ============================================================
#  NORMALIZATION
# ============================================================

"""
    normalizeMatrix(counts::Matrix{Float64}, method::String;
                    meta=nothing, designVar::String="") -> Matrix{Float64}

Normalize a features × samples pseudobulk count matrix.

# Methods
| `method`            | Description                                           | R equivalent          |
|---------------------|-------------------------------------------------------|-----------------------|
| `"none"`            | Raw counts, no normalization                          | `type_norm=NULL`      |
| `"tpm"`             | Transcripts per million                               | `type_norm="TPM"`     |
| `"log2tpm"`         | log2(TPM + 1)                                         | `type_norm="Log2TPM"` |
| `"zscore"`          | Z-score of raw counts per feature across samples      | —                     |
| `"log2tpm_zscore"`  | log2(TPM + 1) then z-score                            | `type_norm="zscore"`  |
| `"log2fc"`          | log2 fold-change relative to row mean                 | `type_norm="Log2FC"`  |
| `"sizefactor"`      | Median-of-ratios size factors + log1p (DESeq2-style)  | —                     |

For exact DESeq2 VST, use `vstNormalizeR` from the
`examples/prepareData/pseudobulkWorkflow.jl` workflow (requires RCall.jl + R with DESeq2).

# Example
```julia
countsNorm = normalizeMatrix(countsMat, "log2tpm")
countsNorm = normalizeMatrix(countsMat, "sizefactor")
```
"""
function normalizeMatrix(counts::Matrix{Float64}, method::String;
                         meta=nothing, designVar::String="")
    if method == "none"
        return counts

    elseif method == "tpm"
        colSums = vec(sum(counts, dims=1))
        colSums[colSums .== 0] .= 1.0
        return 1e6 .* counts ./ colSums'

    elseif method == "log2tpm"
        colSums = vec(sum(counts, dims=1))
        colSums[colSums .== 0] .= 1.0
        tpm = 1e6 .* counts ./ colSums'
        return log2.(tpm .+ 1)

    elseif method == "zscore"
        # Z-score raw counts per feature across samples
        rowMeans = vec(mean(counts, dims=2))
        rowStds  = vec(std(counts,  dims=2))
        rowStds[rowStds .== 0] .= 1.0
        return (counts .- rowMeans) ./ rowStds

    elseif method == "log2tpm_zscore"
        # log2(TPM+1) first, then z-score — equivalent to R's type_norm="zscore"
        log2Tpm  = normalizeMatrix(counts, "log2tpm")
        rowMeans = vec(mean(log2Tpm, dims=2))
        rowStds  = vec(std(log2Tpm,  dims=2))
        rowStds[rowStds .== 0] .= 1.0
        return (log2Tpm .- rowMeans) ./ rowStds

    elseif method == "log2fc"
        colSums = vec(sum(counts, dims=1))
        colSums[colSums .== 0] .= 1.0
        tpm      = 1e6 .* counts ./ colSums' .+ 1
        rowMeans = vec(mean(tpm, dims=2))
        rowMeans[rowMeans .== 0] .= 1.0
        return log2.(tpm ./ rowMeans)

    elseif method == "sizefactor"
        return _sizeFactorNormalize(counts)

    elseif method == "vst"
        error("Exact DESeq2 VST is not available in the package to avoid an R " *
              "dependency. Use `vstNormalizeR` from " *
              "examples/prepareData/pseudobulkWorkflow.jl " *
              "(requires RCall.jl + R with DESeq2 installed).")
    else
        error("Unknown normalization method: \"$method\". " *
              "Choose: none, tpm, log2tpm, zscore, log2tpm_zscore, log2fc, sizefactor")
    end
end


# Internal — called by normalizeMatrix("sizefactor") and directly in workflow
function _sizeFactorNormalize(counts::Matrix{Float64})
    logCounts   = log.(counts .+ 1e-10)
    geomMeans   = vec(mean(logCounts, dims=2))
    logRatios   = logCounts .- geomMeans
    sizeFactors = exp.(vec(median(logRatios, dims=1)))
    sizeFactors[sizeFactors .== 0] .= 1.0
    return log1p.(counts ./ sizeFactors')
end


# ============================================================
#  BATCH CORRECTION
# ============================================================

"""
    removeBatchEffect(X::Matrix{Float64}, batch::Vector;
                      designMat=nothing) -> Matrix{Float64}

Pure-Julia reimplementation of `limma::removeBatchEffect`.

Removes additive batch effects from a features × samples matrix while
preserving biological variation encoded in `designMat`. Fits OLS with
intercept + biological covariates + batch indicators, then subtracts
only the batch component.

Mathematically equivalent to limma. For the exact limma implementation,
use `removeBatchEffectR` from `examples/prepareData/pseudobulkWorkflow.jl`
(requires RCall.jl + R with limma installed).

# Arguments
- `X`          : features × samples normalized matrix
- `batch`      : length-n vector of batch labels (e.g. `["R1","R1","R2",...]`)
- `designMat`  : (optional) n × k biological covariate matrix from
                 `buildDesignMatrix` — preserves these effects during correction

# Example
```julia
bioDesign   = buildDesignMatrix(meta, "celltype")
countsCorrected = removeBatchEffect(countsNorm, meta[!, "replicate"];
                                    designMat=bioDesign)
```
"""
function removeBatchEffect(X::Matrix{Float64}, batch::Vector;
                           designMat::Union{Matrix{Float64}, Nothing}=nothing)
    nSamples  = size(X, 2)
    batches   = sort(unique(batch))
    nBatches  = length(batches)

    batchCols = zeros(Float64, nSamples, nBatches - 1)
    for (j, b) in enumerate(batches[2:end])
        batchCols[:, j] = Float64.(batch .== b)
    end

    intercept  = ones(Float64, nSamples, 1)
    fullDesign = designMat !== nothing ?
                 hcat(intercept, designMat, batchCols) :
                 hcat(intercept, batchCols)
    nProtected = designMat !== nothing ? 1 + size(designMat, 2) : 1

    beta = X * fullDesign * pinv(fullDesign' * fullDesign)
    return X .- beta[:, (nProtected+1):end] * batchCols'
end


"""
    buildDesignMatrix(meta::DataFrame, var::String) -> Matrix{Float64}

One-hot encode a categorical metadata column (dropping the first level)
for use as the biological covariate matrix in `removeBatchEffect`.

# Example
```julia
bioDesign = buildDesignMatrix(meta, "celltype")
```
"""
function buildDesignMatrix(meta::DataFrame, var::String)
    levels = sort(unique(meta[!, var]))
    n      = nrow(meta)
    mat    = zeros(Float64, n, length(levels) - 1)
    for (j, lv) in enumerate(levels[2:end])
        mat[:, j] = Float64.(meta[!, var] .== lv)
    end
    return mat
end


# ============================================================
#  REPLICATE AGGREGATION
# ============================================================

"""
    aggregateReplicates(limaDf::DataFrame,
                        celltypeOrder::Vector{String},
                        replicateOrder::Vector{String};
                        sep::String="_") -> DataFrame

Sum batch-corrected replicate columns into per-celltype aggregate columns.

Column names in `limaDf` must follow `{celltype}{sep}{replicate}`
(e.g. `"gdT17A_R1"`). Returns a DataFrame with `:RowNames` and one
column per cell type present in `celltypeOrder`.

# Example
```julia
limaAgg = aggregateReplicates(limaDf, celltypeOrder, replicateOrder)
```
"""
function aggregateReplicates(limaDf::DataFrame,
                             celltypeOrder::Vector{String},
                             replicateOrder::Vector{String};
                             sep::String="_")
    features   = limaDf[!, :RowNames]
    dataCols   = string.(names(limaDf)[2:end])
    presentCts = [ct for ct in celltypeOrder
                  if any(startswith.(dataCols, ct * sep))]

    aggDf = DataFrame(:RowNames => features)
    for ct in presentCts
        repCols = [c for c in dataCols if startswith(c, ct * sep)]
        if !isempty(repCols)
            aggDf[!, Symbol(ct)] = vec(sum(Matrix(limaDf[:, Symbol.(repCols)]), dims=2))
        end
    end
    return aggDf
end


# vstNormalizeR and removeBatchEffectR are NOT part of this package.
# They are thin wrappers around DESeq2 and limma and belong in user workflow
# scripts, not in a Julia package.
# See examples/prepareData/pseudobulkWorkflow.jl for a complete implementation.
