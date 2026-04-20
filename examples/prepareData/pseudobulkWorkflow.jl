# examples/prepareData/pseudobulkWorkflow.jl
#
# Pseudobulk workflow for InferelatorJL input preparation.
#
# Pipeline
# ────────
#   1. Load AnnData (.h5ad)
#   2. Pseudobulk by replicate-level identity     (e.g. celltype_rep)
#   3. Normalize                                  (DESeq2 VST  OR  package methods)
#   4. Batch correction                           (limma  OR  Julia OLS)
#   5. Save per-replicate outputs
#   6. Pseudobulk by aggregate identity           (e.g. celltype)
#   7. Normalize + aggregate batch-corrected replicates
#   8. Save aggregate outputs
#
# Input  : .h5ad file
#          Convert from Seurat in R:
#            SeuratDisk::SaveH5Seurat(obj, "data.h5seurat")
#            SeuratDisk::Convert("data.h5seurat", dest="h5ad")
#
# Output formats (set via saveFormat below):
#   "tsv"   → tab-delimited text  (human-readable, compatible with R/Python)
#   "arrow" → Apache Arrow        (required input format for InferelatorJL)
#   "both"  → save both formats

using Muon, DataFrames, CSV, Arrow
using RCall          # remove this line if normMethod != "vst" and useRCallBatch = false
using InferelatorJL


# ── RCall helpers (defined here — not part of the InferelatorJL package) ─────

function vstNormalizeR(counts::Matrix{Float64}, meta::DataFrame, designVar::String,
                       blind::Bool)
    @rput counts meta designVar blind
    R"""
    suppressPackageStartupMessages(library(DESeq2))
    countsInt <- round(counts)
    storage.mode(countsInt) <- "integer"
    dds <- DESeqDataSetFromMatrix(
        countData = countsInt,
        colData   = meta,
        design    = as.formula(paste0("~", designVar))
    )
    # vst() requires >= 1000 genes; fall back to the exact method for smaller matrices
    vsd <- tryCatch(
        vst(dds, blind = blind),
        error = function(e) varianceStabilizingTransformation(dds, blind = blind)
    )
    result <- assay(vsd)
    """
    return @rget result
end

function removeBatchEffectR(X::Matrix{Float64}, batch::Vector, designVarVals::Vector)
    @rput X batch designVarVals
    R"""
    suppressPackageStartupMessages(library(limma))
    designMat <- model.matrix(~ designVarVals)
    result    <- limma::removeBatchEffect(X, batch=batch, design=designMat)
    """
    return @rget result
end


# ============================================================
#  USER INPUTS  ← edit this section
# ============================================================

h5adPath    = "/path/to/data.h5ad"
saveDir     = "/path/to/output/PseudobulkData"
assayLabel  = "RNA"             # label used in output filenames ("RNA" or "ATAC")

# Grouping variables (must be columns in adata.obs)
identVar   = "celltype_rep"     # per-replicate grouping
aggVar     = "celltype"         # aggregate-level grouping
batchVar   = "replicate"        # batch variable to remove
designVar  = "celltype"         # biological variable to preserve during batch correction

# Output column order (set to String[] to skip ordering)
celltypeOrder  = ["gdT17A", "gdT17B", "gdT17C", "NaivegdT", "ngdNKT", "gdNKT",
                  "CD8gdT", "ILC3", "ILC1", "ILC2A", "ILC2B", "ILC2C"]
replicateOrder = ["R1", "R2", "R3", "R4"]

# Normalization method — choose one:
#   "none"        raw counts, no normalization
#   "tpm"         transcripts per million
#   "log2tpm"     log2(TPM + 1)
#   "zscore"          z-score of raw counts per feature across samples
#   "log2tpm_zscore"  log2(TPM+1) then z-score  (R's type_norm="zscore")
#   "log2fc"      log2 fold-change relative to row mean
#   "sizefactor"  median-of-ratios size factors + log1p  (pure Julia, no R)
#   "vst"         DESeq2 VST via RCall                   (requires R + DESeq2)
normMethod = "vst"

# VST blind mode:
#   false → use design formula to estimate dispersions (recommended for GRN inference)
#   true  → ignore design, estimate dispersions from all samples
blindBool = false

# Batch correction:
#   true  → limma::removeBatchEffect via RCall  (requires R + limma)
#   false → OLS batch correction in pure Julia  (no R required)
useRCallBatch = true

# Output format:
#   "tsv"   → .tsv files  (human-readable, compatible with R/Python)
#   "arrow" → .arrow files (required input format for InferelatorJL)
#   "both"  → save both
saveFormat = "both"

# ============================================================

mkpath(saveDir)
println("=== Pseudobulk Workflow ===")
println("Norm  : $normMethod + ", useRCallBatch ? "limma (RCall)" : "OLS (Julia)")
println("Input : $h5adPath")
println("Output: $saveDir")


# ── 1. Load ─────────────────────────────────────────────────
println("\n[1/6] Loading AnnData...")
adata = readh5ad(h5adPath)
println("      $(size(adata.X, 1)) cells × $(size(adata.X, 2)) features")


# ── 2. Per-replicate pseudobulk ──────────────────────────────
println("\n[2/6] Pseudobulking by: $identVar")
pbRaw = pseudoBulk(adata, identVar)

allRepCols  = vec([ct * "_" * r for ct in celltypeOrder, r in replicateOrder])
orderedCols = intersect(allRepCols, string.(names(pbRaw)[2:end]))
isempty(orderedCols) && (orderedCols = string.(names(pbRaw)[2:end]))
pbRaw = pbRaw[:, vcat(["RowNames"], orderedCols)]
println("      $(nrow(pbRaw)) features, $(length(orderedCols)) samples")


# ── Per-sample metadata ──────────────────────────────────────
metaWant = intersect([identVar, aggVar, batchVar, designVar], names(adata.obs))
meta     = unique(select(adata.obs, metaWant))
meta     = filter(row -> row[identVar] in orderedCols, meta)
sort!(meta, identVar)

# Derive replicate column from identVar if not already present in obs
# Assumes identVar values follow {celltype}_{replicate} (e.g. "gdT17A_R1" → "R1")
if batchVar ∉ names(meta)
    meta[!, Symbol(batchVar)] = last.(split.(meta[!, identVar], "_"))
end


# ── 3. Normalize ─────────────────────────────────────────────
println("\n[3/6] Normalizing ($normMethod)...")
countsMat  = Matrix(pbRaw[:, 2:end])
countsNorm = normMethod == "vst" ?
             vstNormalizeR(countsMat, meta, designVar, blindBool) :
             normalizeMatrix(countsMat, normMethod)


# ── 4. Batch correction ──────────────────────────────────────
println("\n[4/6] Batch correcting ($batchVar → preserve $designVar)...")
countsLima = if useRCallBatch
    println("      limma::removeBatchEffect (via RCall)")
    removeBatchEffectR(countsNorm, meta[!, batchVar], meta[!, designVar])
else
    println("      OLS batch correction (pure Julia)")
    bioDesign = buildDesignMatrix(meta, designVar)
    removeBatchEffect(countsNorm, meta[!, batchVar]; designMat=bioDesign)
end


# ── 5. Save per-replicate outputs ────────────────────────────
println("\n[5/6] Saving per-replicate outputs...")
prefix   = joinpath(saveDir, "$(assayLabel)_$(designVar)")
features = pbRaw[!, :RowNames]

function _matToDF(mat, colNames, rowNames)
    df = DataFrame(:RowNames => rowNames)
    for (j, cn) in enumerate(colNames)
        df[!, Symbol(cn)] = mat[:, j]
    end
    return df
end

# Save a DataFrame as TSV, Arrow, or both based on saveFormat
# Arrow files rename :RowNames → :Genes to match InferelatorJL input format
function _save(df::DataFrame, basePath::String; format::String=saveFormat)
    if format == "tsv" || format == "both"
        CSV.write(basePath * ".tsv", df; delim='\t')
        println("  Saved → $(basePath).tsv")
    end
    if format == "arrow" || format == "both"
        arrowDf = rename(df, :RowNames => :Genes)
        Arrow.write(basePath * ".arrow", arrowDf)
        println("  Saved → $(basePath).arrow")
    end
end

_save(pbRaw,                                        "$(prefix)_raw")
_save(_matToDF(countsNorm,  orderedCols, features), "$(prefix)_vst")
limaDf = _matToDF(countsLima, orderedCols, features)
_save(limaDf,                                       "$(prefix)_vst_lima")
CSV.write("$(prefix)_metabulk.tsv", meta; delim='\t')   # metadata always TSV
println("  Saved → $(prefix)_metabulk.tsv")


# ── 6. Pseudobulk at aggregate level ─────────────────────────
println("\n[6/6] Aggregating replicates by: $aggVar")
pbAgg      = pseudoBulk(adata, aggVar)
orderedAgg = intersect(celltypeOrder, string.(names(pbAgg)[2:end]))
isempty(orderedAgg) && (orderedAgg = string.(names(pbAgg)[2:end]))
pbAgg = pbAgg[:, vcat(["RowNames"], orderedAgg)]


# ── 7. Normalize aggregated + sum batch-corrected replicates ─
countsAggMat  = Matrix(pbAgg[:, 2:end])
metaAggNorm   = unique(select(adata.obs, intersect([aggVar, designVar], names(adata.obs))))
metaAggNorm   = filter(r -> r[aggVar] in orderedAgg, metaAggNorm)
sort!(metaAggNorm, aggVar)

countsVstAgg = normMethod == "vst" ?
               vstNormalizeR(countsAggMat, metaAggNorm, designVar, blindBool) :
               normalizeMatrix(countsAggMat, normMethod)

limaAgg = aggregateReplicates(limaDf, celltypeOrder, replicateOrder)
limaAgg = limaAgg[:, vcat(["RowNames"], intersect(orderedAgg, string.(names(limaAgg)[2:end])))]

metaAgg = unique(select(adata.obs, [aggVar]))
metaAgg = filter(r -> r[aggVar] in orderedAgg, metaAgg)
sort!(metaAgg, aggVar)


# ── 8. Save aggregate outputs ────────────────────────────────
aggPrefix = joinpath(saveDir, "$(assayLabel)_$(aggVar)")

_save(pbAgg,                                                     "$(aggPrefix)_raw_AggRep")
_save(_matToDF(countsVstAgg, orderedAgg, pbAgg[!, :RowNames]),   "$(aggPrefix)_vst_AggRep")
_save(limaAgg,                                                   "$(aggPrefix)_vst_lima_AggRep")
CSV.write("$(aggPrefix)_metabulk_AggRep.tsv", metaAgg; delim='\t')   # metadata always TSV
println("  Saved → $(aggPrefix)_metabulk_AggRep.tsv")


println("\n✓ Done. All outputs saved to: $saveDir")
