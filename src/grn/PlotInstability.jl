# =============================================================================
#  PlotInstability.jl — Regularization penalty selection diagnostic plots
#  src/grn/PlotInstability.jl
#
#  Generates two figures per pipeline run to validate λ selection:
#
#  Figure 1  instability_diagnostic_<mode>.png   (2 panels side by side)
#  ─────────────────────────────────────────────────────────────────────
#  Left panel  — Model Size vs λ
#    Shows how many regulators appear in the network as λ decreases.
#    At network level : total unique TFs with at least one edge anywhere.
#    At gene level    : # TFs regulating that specific gene.
#    A steep drop as λ increases confirms that regularization is working.
#
#  Right panel — Instability (Lower + Upper Bound) vs λ
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
#  Fig 1 right panel shows the broad warm-start search window.
#  Fig 2 zooms into the fine-pass λ grid within that window and confirms
#  that the chosen λ sits at the minimum of |instability − target|.
#  Together they give a complete picture of the regularisation selection.
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

# Called automatically
Inside `bstartsEstimateInstability` with `mode = :network` after every run.
The figures are saved to the same `outputDir` as `instabOutMat.jld`.

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

    lambdaFine = grnData.lambdaRange       # fine-pass λ grid (x-axis for Fig 2 + model size)
    lambdaWarm = grnData.lambdaRangeWarm   # warm-start λ grid (x-axis for instability bounds)
    totLambdas = length(lambdaFine)

    # ── Select data arrays based on mode ──────────────────────────────────────
    if mode == :network
        # Instability bounds from warm-start pass (one curve per bound)
        instabLb   = grnData.netInstabilitiesLb
        instabUb   = grnData.netInstabilitiesUb

        # Fine-pass instability (one value per λ in lambdaFine)
        instabFine = grnData.netInstabilities

        # Model size: unique TFs with ≥ 1 non-zero edge across all genes at each λ
        modelSize = zeros(Int, totLambdas)
        for lind in 1:totLambdas
            slab = grnData.stabilityMat[lind, :, :]          # totResponses × totPreds
            modelSize[lind] = sum(vec(any(isfinite.(slab) .& (slab .> 0), dims=1)))
        end

        modeLabel   = "Network"
        modelLabel  = "# Unique TF Regulators"
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

        # Model size: # TFs with non-zero edge to this gene at each fine-pass λ
        modelSize = zeros(Int, totLambdas)
        for lind in 1:totLambdas
            row = grnData.stabilityMat[lind, geneIdx, :]
            modelSize[lind] = sum(isfinite.(row) .& (row .> 0))
        end

        modeLabel  = "Gene: $geneName"
        modelLabel = "# TF Regulators"
        figSuffix  = replace(geneName, r"[/\\]" => "_")

    else
        error("mode must be :network or :gene, got :$mode")
    end

    # ── Figure 1: Model size (left) + Instability bounds (right) ─────────────
    fig1, (ax1, ax2) = PyPlot.subplots(1, 2; figsize = (10, 4))
    fig1.suptitle("Regularization Penalty Selection — $modeLabel", fontsize = 13)

    # Left panel — model size vs λ (fine-pass range)
    ax1.plot(lambdaFine, modelSize; color = "steelblue", linewidth = 1.5)
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
        savePath1 = joinpath(outputDir, "instability_diagnostic_$(figSuffix).png")
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
        savePath2 = joinpath(outputDir, "instability_selection_$(figSuffix).png")
        PyPlot.savefig(savePath2; dpi = 150)
        @info "Saved λ selection diagnostic" file = savePath2 chosenλ = chosenλ
    else
        PyPlot.show()
    end
    PyPlot.close(fig2)

    return nothing
end
