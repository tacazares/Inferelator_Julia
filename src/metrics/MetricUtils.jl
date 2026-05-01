"""
    adaptiveStep(xs; target_points=1000, min_step=1e-4, method=:min_gap)

Select a step size for interpolation based on actual values in `xs`.

# Arguments
- `xs`: Array of arrays (per TF) or single array of values.
- `target_points`: Desired number of points if `method=:target_points`.
- `min_step`: Minimum allowable step size.
- `method`: `:min_gap` uses smallest nonzero gap; `:target_points` uses range/target_points.
"""
function adaptiveStep(xs; target_points=1000, min_step=1e-4, method=:min_gap)
    vals = isa(xs[1], AbstractArray) ? reduce(vcat, xs) : xs
    vals = sort(unique(vals))

    if method == :min_gap
        gaps = diff(vals)
        nonzero_gaps = filter(g -> g > 0, gaps)
        step = isempty(nonzero_gaps) ? min_step : max(minimum(nonzero_gaps), min_step)
    elseif method == :target_points
        range_span = maximum(vals) - minimum(vals)
        step = max(range_span / target_points, min_step)
    else
        error("Unknown method: $method")
    end

    return step
end


"""
    dedupNInterpolate(xs, ys, interpPts)

Linearly interpolate `ys` onto `interpPts`, handling duplicates and invalid values.

Filters `NaN`/`Inf`, sorts by `xs`, collapses duplicate `xs` by keeping max `y`,
then interpolates using `Flat()` extrapolation. Returns zeros if fewer than 2 valid points.
"""
function dedupNInterpolate(xs::Vector{Float64}, ys::Vector{Float64}, interpPts::AbstractVector{Float64})
    # Filter out invalid values
    valid_inds = isfinite.(xs) .& isfinite.(ys)
    xs = xs[valid_inds]
    ys = ys[valid_inds]

    # Check if there are enough points to interpolate
    if length(xs) < 2
        return zeros(length(interpPts))
    end

    # Sort xs and corresponding ys
    perm = sortperm(xs)
    xsSorted = xs[perm]
    ysSorted = ys[perm]

    # Collapse duplicates: keep max y for each unique x
    dict = Dict{Float64, Float64}()
    for (x, y) in zip(xsSorted, ysSorted)
        dict[x] = haskey(dict, x) ? max(dict[x], y) : y
    end

    xsUnique = sort(collect(keys(dict)))
    ysUnique = [dict[x] for x in xsUnique]

    # Handle edge case: all ys are identical (LinearInterpolation still works)
    itp = LinearInterpolation(xsUnique, ysUnique, extrapolation_bc=Flat())
    return itp.(interpPts)
end

"""
    validateColumns(cols, required, filename)

Check that all `required` column names are present in `cols`.
Throws an informative error naming the missing columns if any are absent.

usage:
validateColumns(Symbol.(propertynames(gsData)), GS_REQUIRED_COLS, gsFile)
"""
# function validateColumns(cols, required, filename)
#     missing_cols = setdiff(required, cols)
#     if !isempty(missing_cols)
#         error("""
#         Missing columns in $filename: $(join(missing_cols, ", "))
#         Required columns: $(join(required, ", "))
#         Found columns:    $(join(cols, ", "))
#         """)
#     end
# end

function validateColumnCount(data, filename)
    if length(propertynames(data)) < MIN_REQUIRED_COLS
        error("""
        $filename must have at least 3 columns in order: TF, Target, Score
        Found: $(length(propertynames(data))) column(s)
        """)
    end
end


"""
    computeAUPR(recalls, precisions; limit=1.0) -> Float64

Compute the area under the precision-recall curve using the trapezoidal rule.

# Arguments
- `recalls`    : Vector of recall values (x-axis), must be non-decreasing.
- `precisions` : Vector of precision values (y-axis), same length as `recalls`.
- `limit`      : Recall upper bound for partial AUPR (default 1.0 = full AUPR).
                 Only points with recall <= limit are included.

# Returns
Scalar AUPR value. Returns 0.0 if fewer than 2 points fall within the limit.

# Notes
The trapezoidal rule averages adjacent precision values weighted by the recall
width between them. This is more accurate than the step-function approximation
(which uses the left-hand precision for each interval) especially when the
curve is not monotone.

Step-function approximation (kept for reference):
    aupr = 0.0
    prev_recall = 0.0
    for (r, p) in zip(recalls, precisions)
        delta = r - prev_recall
        aupr += delta * p
        prev_recall = r
    end
"""
function computeAUPR(
    recalls::AbstractVector{Float64},
    precisions::AbstractVector{Float64};
    limit::Float64 = 1.0
)::Float64

    if limit < 1.0
        indx       = findall(r -> r <= limit, recalls)
        length(indx) < 2 && return 0.0
        recalls    = recalls[indx]
        precisions = precisions[indx]
    end

    length(recalls) < 2 && return 0.0

    heights = (precisions[2:end] .+ precisions[1:end-1]) ./ 2
    widths  = recalls[2:end] .- recalls[1:end-1]
    return sum(heights .* widths)
end


"""
    computeAUROC(fprs, recalls) -> Float64

Compute the area under the ROC curve using the trapezoidal rule.

# Arguments
- `fprs`    : Vector of false positive rates (x-axis), must be non-decreasing.
- `recalls` : Vector of true positive rates / recall (y-axis), same length as `fprs`.

# Returns
Scalar AUROC value. Returns 0.0 if fewer than 2 points are provided.

# Notes
Step-function approximation (kept for reference):
    aroc = 0.0
    prev_fpr = 0.0
    for (f, r) in zip(fprs, recalls)
        delta = f - prev_fpr
        aroc += delta * r
        prev_fpr = f
    end
"""
function computeAUROC(
    fprs::AbstractVector{Float64},
    recalls::AbstractVector{Float64}
)::Float64

    length(fprs) < 2 && return 0.0

    widths  = fprs[2:end] .- fprs[1:end-1]
    heights = (recalls[2:end] .+ recalls[1:end-1]) ./ 2
    return sum(widths .* heights)
end


"""
    extractMetricsAtLimit(source, limit) -> NamedTuple

Extract four performance metrics at a given recall limit from a saved or
in-memory PR result, without rerunning `computePR`.

# Arguments
- `source` : File path to a saved `.jld` PR result, or an in-memory result
             `Dict` from `computePR` or `loadPRData`.
- `limit`  : Recall limit (e.g. 0.01, 0.05, 0.1). Must be in (0, 1].

# Returns
`NamedTuple` with four fields:

| Field            | Description                                              |
|------------------|----------------------------------------------------------|
| `aupr`           | Partial AUPR over [0, limit] (raw area)                  |
| `auprNormalized` | Partial AUPR / (randPR × limit) — area-based enrichment  |
| `precAtLimit`    | Precision at recall = limit (raw point estimate)         |
| `precNormalized` | precAtLimit / randPR — point-based enrichment            |

All normalized values are `NaN` if `randPR` is zero or unavailable.

# Example
```julia
m = extractMetricsAtLimit(res, 0.05)
m.aupr            # partial AUPR at 5% recall
m.auprNormalized  # area-based enrichment at 5% recall
m.precAtLimit     # precision at 5% recall
m.precNormalized  # point-based enrichment at 5% recall
```
"""
function extractMetricsAtLimit(
    source::Union{String, Dict, OrderedDict},
    limit::Float64
)::NamedTuple

    # Load from file if a path was supplied
    if source isa String
        res = loadPRData(source; mode = :global)
    else
        res = source
    end

    recalls    = res[:recalls]
    precisions = res[:precisions]
    randPR     = get(res, :randPR, NaN)

    # -- Partial AUPR (area-based) --
    aupr           = computeAUPR(recalls, precisions; limit = limit)
    auprNormalized = (randPR > 0 && !isnan(randPR)) ? aupr / (randPR * limit) : NaN

    # -- Precision at limit (point-based) --
    idxAtLimit    = findlast(r -> r <= limit, recalls)
    precAtLimit   = idxAtLimit !== nothing ? precisions[idxAtLimit] : NaN
    precNormalized = (randPR > 0 && !isnan(precAtLimit)) ? precAtLimit / randPR : NaN

    return (
        aupr           = aupr,
        auprNormalized = auprNormalized,
        precAtLimit    = precAtLimit,
        precNormalized = precNormalized
    )
end

"""
    extractMetricsAtLimit(source, limits) -> Dict{Float64, NamedTuple}

Convenience method: compute all four metrics at each limit in `limits`.
Returns a `Dict` mapping each recall limit to its `NamedTuple` of metrics.

# Example
```julia
results = extractMetricsAtLimit(res, [0.01, 0.05, 0.1])
results[0.05].aupr            # partial AUPR at 5% recall
results[0.05].precNormalized  # point-based enrichment at 5% recall
```
"""
function extractMetricsAtLimit(
    source::Union{String, Dict, OrderedDict},
    limits::Vector{Float64}
)::Dict{Float64, NamedTuple}
    # Load once if source is a file path, then reuse the dict for all limits
    res = source isa AbstractString ? loadPRData(source; mode = :global) : source
    return Dict(limit => extractMetricsAtLimit(res, limit) for limit in limits)
end


"""
    saveSummaryTables(tables, gsNames, netNames, dirOut)

Write one TSV summary file per metric table to `dirOut`.

# Arguments
- `tables`   : Vector of `(filename_tag, Dict{gsName => Dict{netName => value}})` pairs.
               Each entry produces one TSV file named `\$(tag)_summary.tsv`.
- `gsNames`  : Ordered vector of gold-standard names (columns in the TSV).
- `netNames` : Ordered vector of network names (rows in the TSV).
- `dirOut`   : Output directory — must already exist.

# Example
```julia
tables = [
    ("aupr",           auprDict),
    ("auprNormalized", auprNormDict),
    ("precAtLimit",    precDict),
    ("precNormalized", precNormDict),
]
saveSummaryTables(tables, gsNames, netNames, dirOutPlot)
```
"""
function saveSummaryTables(
    tables::Vector,
    gsNames::Vector{String},
    netNames::Vector{String},
    dirOut::String
)
    mkpath(dirOut)
    for (tag, dict) in tables
        open(joinpath(dirOut, "$(tag)_summary.tsv"), "w") do f
            println(f, "network\t" * join(gsNames, "\t"))
            for net in netNames
                vals = [string(get(get(dict, gs, OrderedDict()), net, NaN)) for gs in gsNames]
                println(f, net * "\t" * join(vals, "\t"))
            end
        end
    end
    @info "Summary tables saved" dir=dirOut nTables=length(tables)
end