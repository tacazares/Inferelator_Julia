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