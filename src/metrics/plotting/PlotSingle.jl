"""
    plotPRCurve(rec, prec, randPR, saveDir; kwargs...)

Plot and save a Precision-Recall (PR) curve using PyPlot.

# Arguments
- `rec::AbstractVector`        : Recall values.
- `prec::AbstractVector`       : Precision values.
- `randPR::Float64`            : Random baseline precision (horizontal reference line).
- `saveDir::String`            : Directory to save the plot.

# Keyword Arguments
- `xLimitRecall::Float64=1.0`       : Maximum recall shown on the x-axis of diagnostic PR plots.
- `axisTitleSize::Int=16`           : Font size for axis titles.
- `tickLabelSize::Int=14`           : Font size for tick labels.
- `plotTitleSize::Int=18`           : Font size for plot title.
- `baseName=nothing`                : Base filename. If `nothing`, a timestamp is used.
"""
function plotPRCurve(
        rec::AbstractVector, prec::AbstractVector, randPR::Float64, saveDir::String;
        xLimitRecall::Float64 = 1.0, axisTitleSize::Int = 16,
        tickLabelSize::Int = 14, plotTitleSize::Int = 18, baseName = nothing)

    if baseName === nothing
        baseName = Dates.format(now(), "yyyymmdd_HHMMSS")
    end

    if !isempty(rec) && !isempty(prec) && any(prec .> 0)
        PyPlot.figure()
        PyPlot.axhline(y=randPR, linestyle="-.", color="k")
        PyPlot.plot(rec, prec, color="b")
        PyPlot.xlabel("Recall", fontsize=axisTitleSize)
        PyPlot.ylabel("Precision", fontsize=axisTitleSize)
        PyPlot.title("PR Curve: $baseName", fontsize=plotTitleSize)
        PyPlot.xlim(0, xLimitRecall)
        PyPlot.ylim(0, 1)
        PyPlot.grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")
        PyPlot.minorticks_on()
        PyPlot.grid(true, which="minor", linestyle=":", linewidth=0.5, color="lightgray")
        PyPlot.tick_params(axis="both", which="major", labelsize=tickLabelSize)
        PyPlot.tick_params(axis="both", which="minor", labelsize=tickLabelSize - 2)
        PyPlot.savefig(joinpath(saveDir, baseName * "_PR.png"), dpi=600)
        PyPlot.close()
    else
        @warn "Skipping PR plot for $baseName: empty or all-zero precision/recall vectors"
    end
end


"""
    plotROCCurve(fpr, tpr, saveDir; kwargs...)

Plot and save a ROC curve using PyPlot.

# Arguments
- `fpr::AbstractVector` : False positive rate values (x-axis).
- `tpr::AbstractVector` : True positive rate / recall values (y-axis).
- `saveDir::String`     : Directory to save the plot.

# Keyword Arguments
- `axisTitleSize::Int=16`  : Font size for axis titles.
- `tickLabelSize::Int=14`  : Font size for tick labels.
- `plotTitleSize::Int=18`  : Font size for plot title.
- `baseName=nothing`       : Base filename. If `nothing`, a timestamp is used.
"""
function plotROCCurve(
        fpr::AbstractVector, tpr::AbstractVector, saveDir::String;
        axisTitleSize::Int = 16, tickLabelSize::Int = 14,
        plotTitleSize::Int = 18, baseName = nothing)

    if baseName === nothing
        baseName = Dates.format(now(), "yyyymmdd_HHMMSS")
    end

    if !isempty(fpr) && !isempty(tpr) && (maximum(tpr) > 0 || maximum(fpr) > 0)
        PyPlot.figure()
        PyPlot.plot([0, 1], [0, 1], linestyle="--", color="k")
        PyPlot.plot(fpr, tpr, color="b")
        PyPlot.xlabel("FPR", fontsize=axisTitleSize)
        PyPlot.ylabel("TPR", fontsize=axisTitleSize)
        PyPlot.title("ROC Curve: $baseName", fontsize=plotTitleSize)
        PyPlot.xlim(0, 1)
        PyPlot.ylim(0, 1)
        PyPlot.grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")
        PyPlot.minorticks_on()
        PyPlot.grid(true, which="minor", linestyle=":", linewidth=0.5, color="lightgray")
        PyPlot.tick_params(axis="both", which="major", labelsize=tickLabelSize)
        PyPlot.tick_params(axis="both", which="minor", labelsize=tickLabelSize - 2)
        PyPlot.savefig(joinpath(saveDir, baseName * "_ROC.png"), dpi=600)
        PyPlot.close()
    else
        @warn "Skipping ROC plot for $baseName: empty or all-zero TPR/FPR vectors"
    end
end