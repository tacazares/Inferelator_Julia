# Processing/JL/h5adToArrow.jl
#
# Convert an AnnData (.h5ad) expression matrix to Apache Arrow (.arrow)
# for use as input to InferelatorJL.
#
# InferelatorJL expects:
#   - Rows    = genes
#   - Columns = samples / cells
#   - First column named "Genes" containing gene names
#   - Values  = normalized expression (e.g. log-normalized counts)
#
# If starting from a Seurat object, convert to .h5ad first in R:
#   SeuratDisk::SaveH5Seurat(obj, "data.h5seurat")
#   SeuratDisk::Convert("data.h5seurat", dest="h5ad")

using Muon, DataFrames, Arrow


# ============================================================
#  USER INPUTS  ← edit this section
# ============================================================

h5ad_path   = "/path/to/data.h5ad"
output_path = "/path/to/output/expression_logNorm.arrow"

# Which expression matrix to use:
#   nothing  → .X  (the active/default matrix — usually log-normalized)
#   "counts" → raw counts layer
#   "log_norm", "normalized", etc. → any named layer in adata.layers
layer = nothing

# ============================================================

println("Loading: ", h5ad_path)
adata = readh5ad(h5ad_path)

# Extract expression matrix (transpose: cells × features → features × cells)
mat = if layer === nothing
    Matrix(adata.X)'
else
    Matrix(adata.layers[layer])'
end

features = adata.var[!, "name"]         # gene names
cells    = string.(adata.obs_names)     # cell / sample names

println("$(length(features)) genes × $(length(cells)) samples")

# Build DataFrame: Genes column first, then one column per sample
df = DataFrame(mat, Symbol.(cells))
insertcols!(df, 1, :Genes => features)

# Save
mkpath(dirname(output_path))
Arrow.write(output_path, df)
println("Saved → ", output_path)
