# =============================================================================
#  PlotInstability.jl — Regularization penalty selection diagnostic plots
#  src/grn/PlotInstability.jl
#
#  Generates two figures per pipeline run to validate λ selection:
#
#  Figure 1  instability_diagnostic_<mode>.png   (2 panels side by side)
#  ─────────────────────────────────────────────────────────────────────
#  Both panels share the warm-start λ range on the x-axis so model size
#  and instability can be read at the same regularization values.
#
#  Left panel  — Model Size vs λ  (warm-start range)
#    Average # TFs selected per gene across warm-start subsamples.
#    Computed from bstarsWarmStart (theta2 > 0 per gene, averaged).
#    Should decrease monotonically as λ increases — confirms LASSO is
#    working and the warm-start range is wide enough to show variation.
#
#  Right panel — Instability (Lower + Upper Bound) vs λ  (warm-start range)
#    bStARS stability bounds computed during the warm-start pass over the
#    full λ search range [minLambda, maxLambda].
#    The dashed line marks targetInstability (default 0.05).
#    The fine-pass window is the region where Lb ≤ target ≤ Ub.
#
#  Figure 2  instability_selection_<mode>.png    (1 panel)
#  ─────────────────────────────────────────────────────────────────────
#  |netInstabilities − targetInstability| vs λ from the fine estimation pass.
#  The dot marks the chosen λ — the point closest to the target threshold.
#  This is a direct visualisation of what chooseLambda! computes.
#
#  How Figures 1 and 2 relate
#  ─────────────────────────────────────────────────────────────────────
#  Fig 1 shows the broad warm-start search window — model size and
#  instability bounds over the full range, so you can see regularization
#  working. Fig 2 zooms into the fine-pass grid within that window and
#  confirms the chosen λ sits at the minimum of |instability − target|.
#
#  Standalone usage (per-gene diagnostics)
#  ─────────────────────────────────────────────────────────────────────
#  Load the saved GrnData and call plotInstabilityCurves with mode = :gene:
#
#    using JLD2, InferelatorJL
#    grnData = load_object("TFA/instabOutMat.jld")
#    import InferelatorJL: plotInstabilityCurves
#    plotInstabilityCurves(grnData;
#                          mode              = :gene,
#                          geneName          = "IFNG",
#                          targetInstability = 0.05,
#                          outputDir         = "TFA/plots")
# =============================================================================


"""
    plotInstabilityCurves(grnData; mode, geneName, targetInstability, outputDir)

Generate two regularization penalty selection diagnostic figures.

**Figure 1** (`instability_diagnostic_<mode>.png`): two panels showing
model size vs λ (left) and bStARS instability bounds vs λ (right).
The right panel uses the warm-start λ range; the left uses the fine-pass range.

**Figure 2** (`instability_selection_<mode>.png`): |instability − target|
vs λ from the fine estimation pass, with a dot at the selected λ.

# Arguments
- `grnData`            : `GrnData` populated by `bstartsEstimateInstability`
- `mode`               : `:network` (default) — network-level plots using all genes;
                         `:gene` — per-gene plots for one target gene
- `geneName`           : Target gene name string; required when `mode = :gene`
- `targetInstability`  : Instability threshold used during λ selection (default 0.05)
- `outputDir`          : Directory to save both figures; `nothing` displays instead

# On-demand usage
Not called automatically by the pipeline. Load the saved GrnData and call
explicitly when you want to inspect a run's λ selection.

# Standalone (per-gene)
```julia
using JLD2
grnData = load_object("TFA/instabOutMat.jld")
import InferelatorJL: plotInstabilityCurves
plotInstabilityCurves(grnData; mode = :gene, geneName = "IFNG", outputDir = "TFA/plots")
```
"""
function plotInstabilityCurves(
    grnData::GrnData;
    mode::Symbol                      = :network,
    geneName::Union{String, Nothing}  = nothing,
    targetInstability::Float64        = 0.05,
    outputDir::Union{String, Nothing} = nothing
)
    if !isnothing(outputDir)
        mkpath(outputDir)
    end

    lambdaFine = grnData.lambdaRange       # fine-pass λ grid (x-axis for Fig 2)
    lambdaWarm = grnData.lambdaRangeWarm   # warm-start λ grid (x-axis for Fig 1)
    totLambdas = length(lambdaFine)

    # ── Select data arrays based on mode ──────────────────────────────────────
    if mode == :network
        # Instability bounds from warm-start pass (one curve per bound)
        instabLb   = grnData.netInstabilitiesLb
        instabUb   = grnData.netInstabilitiesUb

        # Fine-pass instability (one value per λ in lambdaFine)
        instabFine = grnData.netInstabilities

        # Model size: avg # TFs selected per gene at each warm-start λ
        # (stored during bstarsWarmStart so both panels share the same x-axis)
        modelSize  = grnData.modelSizeWarm

        modeLabel   = "Network"
        modelLabel  = "Avg # TF Regulators per Gene"
        figSuffix   = "network"

    elseif mode == :gene
        isnothing(geneName) && error("geneName must be provided when mode = :gene")
        geneIdx = findfirst(==(geneName), grnData.targGenes)
        isnothing(geneIdx) && error("Gene \"$geneName\" not found in grnData.targGenes")

        # Per-gene instability bounds from warm-start pass
        instabLb   = grnData.instabilitiesLb[geneIdx, :]
        instabUb   = grnData.instabilitiesUb[geneIdx, :]

        # Per-gene fine-pass instability
        instabFine = grnData.geneInstabilities[geneIdx, :]

        # Per-gene model size: # TFs with non-zero stability at each fine-pass λ
        # (fine-pass grid used here because stabilityMat is only on the fine grid)
        modelSizeFine = zeros(Int, totLambdas)
        for lind in 1:totLambdas
            row = grnData.stabilityMat[lind, geneIdx, :]
            modelSizeFine[lind] = sum(isfinite.(row) .& (row .> 0.05))
        end
        modelSize = Float64.(modelSizeFine)

        modeLabel  = "Gene: $geneName"
        modelLabel = "# TF Regulators"
        figSuffix  = replace(geneName, r"[/\\]" => "_")

    else
        error("mode must be :network or :gene, got :$mode")
    end

    # ── Figure 1: Model size (left) + Instability bounds (right) ─────────────
    # Both panels now share lambdaWarm on the x-axis so regularization strength
    # and instability can be read at the same λ values.
    fig1, (ax1, ax2) = PyPlot.subplots(1, 2; figsize = (10, 4))
    fig1.suptitle("Regularization Penalty Selection — $modeLabel", fontsize = 13)

    # Left panel — model size vs λ (warm-start range, same x-axis as right panel)
    lambdaLeft = (mode == :network) ? lambdaWarm : lambdaFine
    ax1.plot(lambdaLeft, modelSize; color = "steelblue", linewidth = 1.5)
    if mode == :network && length(lambdaFine) > 0
        ax1.axvspan(minimum(lambdaFine), maximum(lambdaFine);
                    alpha = 0.15, color = "gray", label = "Fine-pass window")
        ax1.legend(fontsize = 8)
    end
    ax1.set_xscale("log")
    ax1.set_xlabel("λ",          fontsize = 11)
    ax1.set_ylabel(modelLabel,   fontsize = 10)
    ax1.set_title("Model Size vs λ", fontsize = 10)

    # Right panel — instability bounds vs λ (warm-start range)
    ax2.plot(lambdaWarm, instabLb; color = "steelblue", linewidth = 1.5, label = "Lower Bound")
    ax2.plot(lambdaWarm, instabUb; color = "darkorange", linewidth = 1.5, label = "Upper Bound")
    ax2.axhline(y = targetInstability; color = "black", linestyle = "--",
                linewidth = 1.0, label = "Target ($targetInstability)")
    ax2.set_xscale("log")
    ax2.set_xlabel("λ",           fontsize = 11)
    ax2.set_ylabel("Instability", fontsize = 10)
    ax2.set_title("Instability vs λ (warm-start range)", fontsize = 10)
    ax2.legend(fontsize = 8)

    PyPlot.tight_layout()

    if !isnothing(outputDir)
        savePath1 = joinpath(outputDir, "warmstart_bounds_$(figSuffix).png")
        PyPlot.savefig(savePath1; dpi = 150)
        @info "Saved instability diagnostic" file = savePath1
    else
        PyPlot.show()
    end
    PyPlot.close(fig1)

    # ── Figure 2: |instability − target| with chosen λ dot ───────────────────
    devs    = abs.(instabFine .- targetInstability)
    minIdx  = argmin(devs)
    chosenλ = lambdaFine[minIdx]

    fig2, ax = PyPlot.subplots(1, 1; figsize = (5, 4))
    fig2.suptitle("λ Selection — $modeLabel", fontsize = 13)

    ax.plot(lambdaFine, devs; color = "steelblue", linewidth = 1.5)
    ax.scatter([chosenλ], [devs[minIdx]]; color = "steelblue", s = 80, zorder = 5)
    ax.set_xscale("log")
    ax.set_xlabel("λ",                       fontsize = 11)
    ax.set_ylabel("|Instability − Target|",  fontsize = 10)
    ax.set_title("|Instability − Target| vs λ", fontsize = 10)

    PyPlot.tight_layout()

    if !isnothing(outputDir)
        savePath2 = joinpath(outputDir, "lambda_selection_$(figSuffix).png")
        PyPlot.savefig(savePath2; dpi = 150)
        @info "Saved λ selection diagnostic" file = savePath2 chosenλ = chosenλ
    else
        PyPlot.show()
    end
    PyPlot.close(fig2)
    return nothing
end
