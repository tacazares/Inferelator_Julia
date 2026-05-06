using CSV, DataFrames, JLD2, Statistics, StatsBase, PyPlot, Dates, OrderedCollections

# Set global font defaults (matches plotEdgeConfidence.jl)
PyPlot.matplotlib.rcParams["font.family"]    = "Nimbus Sans"
PyPlot.matplotlib.rcParams["axes.titlesize"] = 12
PyPlot.matplotlib.rcParams["axes.labelsize"] = 10
PyPlot.matplotlib.rcParams["xtick.labelsize"] = 8
PyPlot.matplotlib.rcParams["ytick.labelsize"] = 8
PyPlot.matplotlib.rcParams["legend.fontsize"] = 9


# =============================================================================
#  validateBEBICLambdas.jl — Compare per-gene EBIC-optimal lambda distributions
#                            across bEBIC runs (auto vs fixed-step grids)
#
#  Confirms that user-supplied lambda grids (40 or 100 steps) resolve the same
#  lambda region as GLMNet's automatic path. Overlapping ECDFs = equivalent grids.
#
#  Three data sources (set per run via `source`):
#
#    :tsv        — reads bebic_lambda_summary.tsv (gene | median_lambda | std_lambda)
#                  Fast. Only available for runs where bebicFinalFit! was called
#                  with outputDir set.
#
#    :jld        — tries bebicOutMat.jld first (saves full grnData: gene names +
#                  ebicLambdas + lambdaSS → can compute std), then falls back to
#                  grnOutMat.jld (saves buildGrn: median lambda only, no gene
#                  names, no std).
#
#    :auto       — tries TSV first, then bebicOutMat.jld, then grnOutMat.jld.
#                  Picks the richest available source automatically.
#
#  Outputs (saved to dirOut):
#    lambdaECDF_medianLambda_<subfolder>.pdf   — ECDF of per-gene median lambda
#    lambdaDensity_stdLambda_<subfolder>.pdf   — density of lambda std (where available)
#    lambdaSummary_stats_<subfolder>.tsv       — numeric summary table
# =============================================================================


# =============================================================================
# USER INPUTS
# =============================================================================

# Named dict: label => (dir = path to TFA or TFmRNA subfolder, source = :auto/:tsv/:jld)
runDefs = OrderedDict(
    "bEBIC_autoLambda"    => (dir    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max/TFmRNA",
                               source = :jld),   # has bebicOutMat.jld → full lambdaSS
    "bEBIC_40lambda"      => (dir    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_lambdaProvided/TFmRNA",
                               source = :jld),   # only grnOutMat.jld → median only
    "bEBIC_100lambda"     => (dir    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_100lambdaProvided/TFmRNA",
                               source = :tsv),   # has bebic_lambda_summary.tsv
    "bEBIC_noZ_40lambda"  => (dir    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_noZscoreTFA_noZscoreLASSO_lambdaProvided/TFmRNA",
                               source = :jld),   # only grnOutMat.jld → median only
    "bEBIC_noZ_100lambda" => (dir    = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/bEBIC_gamma0p5_80totSS_subsamplePCT68_0p5_20tfsPerGene_max_noZscoreTFA_noZscoreLASSO_100lambdaProvided/TFmRNA",
                               source = :auto),  # has both TSV and bebicOutMat.jld
)

# Subfolder label used in output filenames
subfolderLabel = "TFmRNA"

# Output directory
dirOut = "/Users/owop7y/Desktop/InferelatorJL/outputs/Bulk/netEval/lambdaValidation"
mkpath(dirOut)


# =============================================================================
# LOADERS
# =============================================================================

"""
Load lambda data from bebic_lambda_summary.tsv.
Returns (genes, medianLambda, stdLambda) — all fields available.
"""
function _loadTSV(dir::String)
    path = joinpath(dir, "bebic_lambda_summary.tsv")
    if !isfile(path)
        @warn "TSV not found: $path"
        return nothing
    end
    df = CSV.read(path, DataFrame; delim='\t')
    return (genes       = df.gene,
            medianLambda = Float64.(df.median_lambda),
            stdLambda    = Float64.(df.std_lambda),
            source       = "TSV")
end

"""
Load lambda data from bebicOutMat.jld (saves full grnData).
Returns (genes, medianLambda, stdLambda) with std computed from lambdaSS.
"""
function _loadBebicJLD(dir::String)
    path = joinpath(dir, "bebicOutMat.jld")
    if !isfile(path)
        return nothing
    end
    grnData = load_object(path)
    if isempty(grnData.ebicLambdas)
        @warn "bebicOutMat.jld at $dir has empty ebicLambdas"
        return nothing
    end
    stdLambda = if !isempty(grnData.lambdaSS) && size(grnData.lambdaSS, 2) > 1
        [std(grnData.lambdaSS[i, :]) for i in 1:size(grnData.lambdaSS, 1)]
    else
        fill(NaN, length(grnData.ebicLambdas))
    end
    return (genes        = grnData.targGenes,
            medianLambda = Float64.(grnData.ebicLambdas),
            stdLambda    = stdLambda,
            source       = "bebicOutMat.jld")
end

"""
Load lambda data from grnOutMat.jld (saves buildGrn — median lambda only).
No gene names, no std. Used as last-resort fallback.
"""
function _loadGrnJLD(dir::String)
    path = joinpath(dir, "grnOutMat.jld")
    if !isfile(path)
        @warn "grnOutMat.jld not found: $path"
        return nothing
    end
    buildGrn = load_object(path)
    lambdaVec = vec(buildGrn.lambda)
    if isempty(lambdaVec) || all(lambdaVec .== 0)
        @warn "grnOutMat.jld at $dir has empty/zero lambda — bStARS run? Skipping."
        return nothing
    end
    n = length(lambdaVec)
    @info "  Loaded grnOutMat.jld ($n genes, median-only, no std)"
    return (genes        = ["gene_$i" for i in 1:n],   # no gene names in buildGrn
            medianLambda = lambdaVec,
            stdLambda    = fill(NaN, n),
            source       = "grnOutMat.jld (median only)")
end

"""
Load lambda data for one run according to the requested source strategy.
  :tsv  — TSV only
  :jld  — bebicOutMat.jld → grnOutMat.jld
  :auto — TSV → bebicOutMat.jld → grnOutMat.jld
"""
function loadLambdaData(dir::String, source::Symbol)
    if source == :tsv
        return _loadTSV(dir)

    elseif source == :jld
        result = _loadBebicJLD(dir)
        if !isnothing(result); return result; end
        @info "bebicOutMat.jld not found in $dir — falling back to grnOutMat.jld (median only)"
        return _loadGrnJLD(dir)

    elseif source == :auto
        result = _loadTSV(dir)
        if !isnothing(result); return result; end
        result = _loadBebicJLD(dir)
        if !isnothing(result); return result; end
        @info "Neither TSV nor bebicOutMat.jld in $dir — falling back to grnOutMat.jld (median only)"
        return _loadGrnJLD(dir)

    else
        error("Unknown source :$source. Use :tsv, :jld, or :auto.")
    end
end


# =============================================================================
# LOAD ALL RUNS
# =============================================================================

runData = OrderedDict{String, Any}()

for (nm, def) in runDefs
    @info "Loading: $nm  [source=:$(def.source)]  $(def.dir)"
    result = loadLambdaData(def.dir, def.source)
    if isnothing(result)
        @warn "  Skipping $nm — no data available."
        continue
    end
    @info "  → $(length(result.medianLambda)) genes  |  source: $(result.source)"
    runData[nm] = result
end

if isempty(runData)
    error("No runs loaded. Check paths and source settings.")
end

# Runs that have valid std (not all-NaN)
runsWithStd = filter(nm -> !all(isnan.(runData[nm].stdLambda)), keys(runData))


# =============================================================================
# COLOUR MAP  (consistent across all plots)
# =============================================================================

using PyPlot: matplotlib
cmap   = matplotlib.cm.get_cmap("tab10")
colors = Dict(nm => cmap((i-1) / max(length(runData)-1, 1)) for (i, nm) in enumerate(keys(runData)))


# =============================================================================
# PLOT 1 — ECDF of per-gene median lambda (log scale)
#
# Overlapping curves → grids resolve the same lambda region → equivalent.
# A horizontally shifted curve → grid boundaries were too narrow.
# =============================================================================

fig1, ax1 = subplots(figsize=(6, 4), layout="constrained")

for nm in keys(runData)
    vals    = runData[nm].medianLambda
    vals    = filter(x -> x > 0 && isfinite(x), vals)
    sorted  = sort(vals)
    n       = length(sorted)
    ecdf_y  = (1:n) ./ n
    ax1.plot(sorted, ecdf_y; label=nm, color=colors[nm], linewidth=1.2)
end

ax1.set_xscale("log")
ax1.set_xlabel("Per-gene median lambda (log scale)")
ax1.set_ylabel("Cumulative proportion of genes")
ax1.set_title("bEBIC lambda ECDF — $subfolderLabel")
ax1.set_ylim(0, 1)
ax1.set_yticks(0:0.2:1)
ax1.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))

# Compute x-axis limits that span full decades so at least 2 decade labels
# appear even when all runs' lambdas fall within a single decade.
_all_ecdf_vals = vcat([filter(x -> x > 0 && isfinite(x), runData[nm].medianLambda)
                       for nm in keys(runData)]...)
_xlo = 10.0 ^ floor(log10(minimum(_all_ecdf_vals)))
_xhi = 10.0 ^ ceil(log10(maximum(_all_ecdf_vals)))
# If data spans < 2 decades, expand by one decade on each side so multiple
# decade labels are always visible.
if log10(_xhi) - log10(_xlo) < 2
    _xlo /= 10.0
    _xhi *= 10.0
end
ax1.set_xlim(_xlo, _xhi)

ax1.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=15))
ax1.xaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation(base=10, labelOnlyBase=true))
ax1.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10, subs=collect(2:9) ./ 10, numticks=100))
ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax1.grid(true, which="major", linestyle="-",  linewidth=0.5, color="lightgray")
ax1.grid(true, which="minor", linestyle=":",  linewidth=0.2, color="#e8e8e8")
ax1.legend(frameon=false, loc="lower right")

out1 = joinpath(dirOut, "lambdaECDF_medianLambda_$(subfolderLabel).pdf")
PyPlot.savefig(out1, dpi=600)
PyPlot.close(fig1)
println("Saved: $out1")


# =============================================================================
# PLOT 2 — Density of per-gene lambda std (log scale)
#
# High std → EBIC-optimal lambda varies widely across subsamples for that gene.
# Only plotted for runs where std is available (bebicOutMat.jld or TSV).
# Expected: z-scored runs have lower std (lambda scale normalised across genes).
# =============================================================================

if !isempty(runsWithStd)
    fig2, ax2 = subplots(figsize=(6, 4), layout="constrained")

    for nm in runsWithStd
        vals = runData[nm].stdLambda
        vals = filter(x -> x > 0 && isfinite(x), vals)   # drop 0-std and NaN
        if isempty(vals)
            @warn "  No valid std values for $nm — skipping density"
            continue
        end
        # Kernel density estimate on log scale
        logVals = log10.(vals)
        kde     = fit(UnivariateKDE, logVals)
        xs      = 10 .^ kde.x
        ax2.plot(xs, kde.density ./ (xs .* log(10));  # Jacobian for log scale
                 label=nm, color=colors[nm], linewidth=1.2)
        ax2.fill_between(xs, kde.density ./ (xs .* log(10)); alpha=0.12, color=colors[nm])
    end

    ax2.set_xscale("log")
    ax2.set_xlabel("Per-gene lambda std across subsamples (log scale)")
    ax2.set_ylabel("Density")
    ax2.set_title("bEBIC lambda variability — $subfolderLabel")
    ax2.minorticks_on()
    ax2.grid(true, which="major", linestyle="-",  linewidth=0.5, color="lightgray")
    ax2.grid(true, which="minor", linestyle=":",  linewidth=0.2, color="#e8e8e8")
    ax2.legend(frameon=false)

    out2 = joinpath(dirOut, "lambdaDensity_stdLambda_$(subfolderLabel).pdf")
    PyPlot.savefig(out2, dpi=600)
    PyPlot.close(fig2)
    println("Saved: $out2")

    if length(runsWithStd) < length(runData)
        missing_std    = setdiff(collect(keys(runData)), collect(runsWithStd))
        missing_str    = join(missing_std, ", ")
        @warn "Std not available for: $missing_str\n" *
              "  → Run bebicFinalFit!(grnData, buildGrn; outputDir=dir) to save bebicOutMat.jld"
    end
else
    @warn "No runs have std data available — density plot skipped."
end


# =============================================================================
# SUMMARY TABLE
# =============================================================================

println("\n", "="^70)
println("Lambda summary — $subfolderLabel")
println("="^70)
println(rpad("Run", 25), rpad("n", 7), rpad("min", 10), rpad("Q25", 10),
        rpad("median", 10), rpad("Q75", 10), rpad("max", 10), "std_median")
println("-"^70)

summaryRows = []
for nm in keys(runData)
    vals = filter(isfinite, runData[nm].medianLambda)
    std_med = let s = filter(x -> isfinite(x) && x > 0, runData[nm].stdLambda)
        isempty(s) ? NaN : median(s)
    end
    row = (run       = nm,
           n_genes   = length(vals),
           λ_min     = minimum(vals),
           λ_q25     = quantile(vals, 0.25),
           λ_median  = median(vals),
           λ_q75     = quantile(vals, 0.75),
           λ_max     = maximum(vals),
           std_median = std_med)
    push!(summaryRows, row)
    println(rpad(nm, 25),
            rpad(row.n_genes, 7),
            rpad(round(row.λ_min,    sigdigits=3), 10),
            rpad(round(row.λ_q25,    sigdigits=3), 10),
            rpad(round(row.λ_median, sigdigits=3), 10),
            rpad(round(row.λ_q75,    sigdigits=3), 10),
            rpad(round(row.λ_max,    sigdigits=3), 10),
            isnan(row.std_median) ? "N/A" : round(row.std_median, sigdigits=3))
end
println("="^70)

# Save summary TSV
summaryDf = DataFrame(summaryRows)
outTsv    = joinpath(dirOut, "lambdaSummary_stats_$(subfolderLabel).tsv")
CSV.write(outTsv, summaryDf; delim='\t')
println("Saved: $outTsv")
