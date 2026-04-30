# =============================================================================
#  plotPR.jl — Evaluate GRN predictions against gold standards and plot PR curves
#
#  What:
#    Evaluates one or more inferred GRNs against gold-standard interaction sets
#    by computing precision-recall (PR) and ROC metrics, then generating
#    publication-quality PR curve plots. Supports multiple networks and
#    multiple gold standards in a single run. Optionally generates per-TF
#    PR curves and AUPR bar plots.
#
#  Required inputs:
#    outNetFiles    — inferred GRN file(s) as legend label → file path dict
#    gsParam        — gold-standard file(s) as name → file path dict
#    prTargGeneFile — target gene list used to restrict evaluation universe
#                     (set to "" to use all genes in the network)
#    gsRegsFile     — regulator list to restrict evaluation to shared TFs
#                     (set to "" to use all regulators)
#
#  Expected outputs (written relative to each network file's directory):
#    PR_noPotRegs/<gsName>/          — PR data files (if gsRegsFile = "")
#    PR_withPotRegs/<gsName>/        — PR data files (if gsRegsFile is set)
#    dirOutPlot/<figBaseName>_*.png  — PR curve plots and optional AUPR bar plots
#
#  Usage:
#    julia examples/plotPR.jl
#    or configure the USER CONFIG section and run step-by-step in the REPL
#
#  Installation:
#    pkg> dev /path/to/InferelatorJL      # local development
#    pkg> add "https://github.com/org/InferelatorJL.jl"   # published release
#    Tip: load Revise before this file to pick up source edits without restarting:
#         using Revise; using InferelatorJL
# =============================================================================

using Revise
using InferelatorJL
import InferelatorJL: computePR, plotPRCurves, plotAUPR, loadPRData,
                      extractMetricsAtLimit, saveSummaryTables

using OrderedCollections

# =========================================================
#  USER CONFIG — edit this section
# =========================================================

# Output directory for plots
dirOutPlot = "/path/to/plots"

# Base name for saved figures (set to "" to use gold-standard name only)
figBaseName = "myRun"

# Network files to compare: legend label => file path
# Add one entry per network you want to compare
outNetFiles = OrderedDict(
    "TFA"      => "/path/to/output/networkLambda.../TFA/edges_subset.tsv",
    "TFmRNA"   => "/path/to/output/networkLambda.../TFmRNA/edges_subset.tsv",
    "Combined" => "/path/to/output/networkLambda.../Combined/combined_max.tsv",
)

# Gold-standard files: name => file path
# Add one entry per gold standard you want to evaluate against
gsParam = OrderedDict(
    "gsName" => "/path/to/goldStandard.tsv",
)

# Evaluation inputs
prTargGeneFile = "/path/to/target_genes.txt"   # set to "" to use all genes
gsRegsFile     = "/path/to/potential_regs.txt" # set to "" to use all TFs
breakTies      = true
auprLimit     = 0.1              # partial AUPR limit passed to computePR
summaryLimits = [0.01, 0.05, 0.1] # recall limits for summary tables — change freely, never reruns computePR

# Plot parameters
lineTypes    = []   # e.g. ["-", "--", "-."] — one per dataset; [] uses defaults
lineWidths   = []   # per dataset; [] uses defaults
lineColors   = []   # per dataset; [] uses defaults
xLimitRecall = 0.1
yStepSize    = [0.1, 0.1]  # one per gold standard
yScaleType   = "linear"
yZoomPR      = [[0.3, 0.9], [], []]  # one per gold standard, or []
heightRatios = [0.5, 3.0]  # height ratios for broken y-axis panels
isInside     = false   # legend inside the plot
plotAUPRflag = false   # set to true to also generate AUPR bar plots
combinePlot  = true    # generate a combined PR curve per network/GS pair
doPerTF      = true    # compute per-TF PR metrics
tfList       = []      # list of TF names for per-TF curves; [] skips Section 3

# =========================================================
#  EXECUTION — no edits needed below this line
# =========================================================

mkpath(dirOutPlot)

# --- Helper: resolve per-GS plot parameters ---
function getPlotParams(i, gsName; figBaseName, yZoomPR, yStepSize)
    saveNamePR       = isempty(figBaseName) ? "$(gsName)" : "$(figBaseName)_$(gsName)"
    currentYzoomPR   = (length(yZoomPR) >= i && !isempty(yZoomPR[i])) ? yZoomPR[i] : Float64[]
    currentYstepSize = (length(yStepSize) >= i && !isempty(yStepSize[i])) ? yStepSize[i] : nothing
    return saveNamePR, currentYzoomPR, currentYstepSize
end

# ----- 1A. Compute PR/ROC metrics and save results ------------------------------
# Run once — results are saved to .jld files in each network's PR directory.
# Comment out this section on subsequent runs to skip recomputation.
@info "---- 1. Calculating Performance Metrics for the Networks -----"

parts = String[]
isempty(gsRegsFile)     || push!(parts, "regs")
isempty(prTargGeneFile) || push!(parts, "targs")
prSuffix = isempty(parts) ? "" : "_" * join(parts, "_")

# prFilesByGS stores saved .jld file paths — required by plotPRCurves and plotAUPR.
# To skip Section 1 on re-runs, populate manually:
#   prFilesByGS = OrderedDict("gsName" => OrderedDict("TFA" => "/path/to/.jld", ...))
prFilesByGS = OrderedDict{String, OrderedDict{String, String}}()

for (legendLabel, outNetFile) in outNetFiles
    @info "Processing network" network=legendLabel file=outNetFile
    filepath = dirname(outNetFile)

    for (gsName, gsFile) in gsParam
        dirPR = joinpath(filepath, "PR" * prSuffix, gsName)
        mkpath(dirPR)
        @info "Using GS" gs=gsName saveDir=dirPR

        res = computePR(gsFile, outNetFile;
                        gsRegsFile       = gsRegsFile,
                        targGeneFile     = prTargGeneFile,
                        breakTies        = breakTies,
                        partialAUPRlimit = auprLimit,
                        doPerTF          = doPerTF,
                        saveDir          = dirPR)

        savedFile = res[:savedFile]
        if savedFile === nothing
            @warn "computePR returned no saved file — skipping" network=legendLabel gs=gsName
            continue
        end

        !haskey(prFilesByGS, gsName) && (prFilesByGS[gsName] = OrderedDict{String, String}())
        prFilesByGS[gsName][legendLabel] = savedFile
    end
end

# ----- 1B. Extract summary metrics at multiple recall limits ------------------
# Fast: loads saved .jld files — re-run freely without rerunning Section 1.
# Adjust summaryLimits in USER CONFIG and re-run only this section.
@info "---- 1.5. Extracting Summary Metrics -----"

netNames = collect(keys(outNetFiles))
gsNames  = collect(keys(gsParam))

auprFull        = OrderedDict{String, OrderedDict{String, Float64}}()
auprAtLimit     = Dict(l => OrderedDict{String, OrderedDict{String, Float64}}() for l in summaryLimits)
auprNormAtLimit = Dict(l => OrderedDict{String, OrderedDict{String, Float64}}() for l in summaryLimits)
precAtLimit     = Dict(l => OrderedDict{String, OrderedDict{String, Float64}}() for l in summaryLimits)
precNormAtLimit = Dict(l => OrderedDict{String, OrderedDict{String, Float64}}() for l in summaryLimits)

for (gsName, netDict) in prFilesByGS
    auprFull[gsName] = OrderedDict{String, Float64}()
    for l in summaryLimits
        auprAtLimit[l][gsName]     = OrderedDict{String, Float64}()
        auprNormAtLimit[l][gsName] = OrderedDict{String, Float64}()
        precAtLimit[l][gsName]     = OrderedDict{String, Float64}()
        precNormAtLimit[l][gsName] = OrderedDict{String, Float64}()
    end

    for (legendLabel, source) in netDict
        # include 1.0 to get full AUPR in the same pass — single load regardless of source type
        allM = extractMetricsAtLimit(source, vcat(summaryLimits, [1.0]))
        auprFull[gsName][legendLabel] = allM[1.0].aupr

        for l in summaryLimits
            m = allM[l]
            auprAtLimit[l][gsName][legendLabel]     = m.aupr
            auprNormAtLimit[l][gsName][legendLabel] = m.auprNormalized
            precAtLimit[l][gsName][legendLabel]     = m.precAtLimit
            precNormAtLimit[l][gsName][legendLabel] = m.precNormalized
        end
    end
end

tables = [("aupr_full", auprFull)]
for l in summaryLimits
    push!(tables, ("aupr_partial_$(l)",    auprAtLimit[l]))
    push!(tables, ("aupr_normalized_$(l)", auprNormAtLimit[l]))
    push!(tables, ("prec_at_limit_$(l)",   precAtLimit[l]))
    push!(tables, ("prec_normalized_$(l)", precNormAtLimit[l]))
end

saveSummaryTables(tables, gsNames, netNames, dirOutPlot)

# Example filter
# keepKeys = ["eq1.gMax", "ge05lt1.gMax"]   # replace with the keys you want
# subset = OrderedDict(
#     gs => OrderedDict(k => v for (k, v) in nets if k in keepKeys)
#     for (gs, nets) in aa
# )
# prFilesByGS = subset

# ----- 2. Global PR curves --------------------------------------------
@info "---- 2. Generating Global PR Curves ----"
if combinePlot
    for (i, (gsName, listFilePR)) in enumerate(prFilesByGS)
        @info "Plotting PR curves" gs=gsName

        saveNamePR, currentYzoomPR, currentYstepSize = getPlotParams(i, gsName;
                                                                      figBaseName = figBaseName,
                                                                      yZoomPR     = yZoomPR,
                                                                      yStepSize   = yStepSize)

        plotPRCurves(listFilePR, dirOutPlot, saveNamePR;
                     xLimitRecall = xLimitRecall,
                     yZoomPR      = currentYzoomPR,
                     yStepSize    = currentYstepSize,
                     yScale       = yScaleType,
                     isInside     = isInside,
                     lineColors   = lineColors,
                     lineTypes    = lineTypes,
                     lineWidths   = lineWidths,
                     heightRatios = heightRatios,
                     mode         = :global)

        if plotAUPRflag
            singleGS     = OrderedDict(gsName => listFilePR)
            saveNameAUPR = isempty(figBaseName) ? "$(gsName)" : "$(figBaseName)_$(gsName)"

            for (figSize, saveLegend) in [((5, 4), true), ((1.5, 1.5), false)]
                plotAUPR(singleGS, dirOutPlot;
                         saveName       = saveNameAUPR,
                         metricType     = "partial",
                         figSize        = figSize,
                         axisTitleSize  = 9,
                         tickLabelSize  = 7,
                         legendFontSize = 9,
                         tickRotation   = 45,
                         plotType       = "bar",
                         saveLegend     = saveLegend)
            end
        end
        @info "Plots completed" gs=gsName
    end
end

# ----- 3. Per-TF PR curves --------------------------------------------
if !isempty(tfList)
    @info "----- 3. Generating Per-TF PR Curves -----"
    for (i, (gsName, resultsDict)) in enumerate(prFilesByGS)
        @info "Plotting per-TF PR curves" gs=gsName

        saveNamePR, currentYzoomPR, currentYstepSize = getPlotParams(i, gsName;
                                                                      figBaseName = figBaseName,
                                                                      yZoomPR     = yZoomPR,
                                                                      yStepSize   = yStepSize)
        saveNamePR = "perTF_$(saveNamePR)"
        tfListPR   = OrderedDict()

        resCache = Dict{String, Any}()
        for (runName, source) in resultsDict
            resCache[runName] = loadPRData(source; mode = :perTF)
        end

        for (runName, res) in resCache
            res === nothing && continue
            tfIndex = Dict(tf => j for (j, tf) in enumerate(res[:gsRegs]))
            for tf in tfList
                idx = get(tfIndex, tf, nothing)
                idx === nothing && continue

                label =
                    length(tfList) == 1 && length(resultsDict) > 1 ? runName :
                    length(resultsDict) == 1 && length(tfList) > 1 ? tf :
                    "$runName - $tf"

                tfListPR[label] = Dict(
                    :precisions => res[:precisions][idx],
                    :recalls    => res[:recalls][idx],
                    :randPR     => res[:randPR][idx]
                )
            end
        end

        plotPRCurves(tfListPR, dirOutPlot, saveNamePR;
                     xLimitRecall = xLimitRecall,
                     yZoomPR      = currentYzoomPR,
                     yStepSize    = currentYstepSize,
                     yScale       = yScaleType,
                     isInside     = isInside,
                     lineColors   = lineColors,
                     lineTypes    = lineTypes,
                     lineWidths   = lineWidths)
    end
end

@info "Completed — plots generated for all gold standards"
