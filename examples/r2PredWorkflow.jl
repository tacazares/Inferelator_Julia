# examples/r2PredWorkflow.jl
#
# Out-of-sample R² prediction workflow for model-size selection.
#
# Uses leave-out cross-validation to select the optimal number of TFs per
# gene (model size) for network inference with mLASSO-StARS.
#
# Pipeline
# ────────
#   1. Load expression data, prior, estimate TFA  (standard pipeline steps 1–3)
#   2. For each leave-out set × predictor mode (TFA / TF mRNA):
#        a. Build GRN holding out the leave-out samples from subsampling
#        b. Load saved GrnData (written by bstartsEstimateInstability)
#        c. Compute R²_pred across model sizes 1:meanEdgesPerGene
#        d. Save r2_pred.pdf, model_size.pdf, r2pred.jld2
#
# References
# ──────────
#   Miraldi et al. (2019) PLOS Comput. Biol.  — mLASSO-StARS
#   Liu, Roeder, Wasserman (2010)             — StARS
#   Muller, Kurtz, Bonneau (2016) arXiv       — bStARS
#
# Output
# ──────
#   For each leave-out × mode combination, outputs are written to:
#     <outputDir>/<leaveOutLabel>/<mode>/
#       edges.tsv, edges_subset.tsv    ← ranked network
#       instabOutMat.jld               ← serialised GrnData (loaded for R²)
#       r2_pred.pdf                    ← R²_pred / R²_fit vs instability cutoff
#       model_size.pdf                 ← R²_pred vs # TFs per gene

using InferelatorJL
using JLD2    # for load_object (reloads GrnData written by buildNetwork)

# ── USER INPUTS ──────────────────────────────────────────────────────────────

normGeneExprFile   = "/path/to/expression.arrow"   # genes × samples (Arrow or TSV)
targGeneFile       = "/path/to/target_genes.txt"   # one gene per line
potRegFile         = "/path/to/pot_regs.txt"       # one TF per line
priorFile          = "/path/to/prior.tsv"
priorFilePenalties = ["/path/to/prior.tsv"]        # same as prior or separate penalty file

# Leave-out sets: each entry is [path_to_sample_list, label]
# The sample list is a plain-text file with one sample name per line.
loInfo = [
    ["/path/to/leaveout_set1.tsv", "ConditionA"],
    ["/path/to/leaveout_set2.tsv", "ConditionB"],
]

outputDir    = "/path/to/outputs/R2Pred"
tfaOptions   = [true, false]   # true = TFA predictors, false = TF mRNA predictors
modeLabels   = ["TFA", "TFmRNA"]

# ── Parameters ───────────────────────────────────────────────────────────────

modelSelection   = "bStARS"    # "bStARS" or "bEBIC" — both supported for R² evaluation
                               # "EBIC" is not supported (no leave-out split, no selection frequencies)

# --- Subsampling and stability
totSS            = 80          # total bootstrap subsamples for stability estimation
bstarsTotSS      = 5           # subsamples for warm-start lambda range search (coarser, faster)
subsampleFrac    = 0.63        # fraction of samples per subsample (0.63–0.68 typical)
targetInstability = 0.05       # instability threshold for bStARS lambda selection (0.05 = 5%)
instabilityLevel = "Network"   # "Network": one shared lambda for all genes
                               # "Gene"   : per-gene lambda — slower, more flexible

# --- Lambda grid
lambdaMin        = 0.01        # LASSO lambda search range — lower bound
lambdaMax        = 1.0         # LASSO lambda search range — upper bound
totLambdasBstars = 20          # lambdas tested in warm-start phase
totLambdas       = 40          # lambdas tested in full stability estimation

# --- Network structure
meanEdgesPerGene = 20          # average TF regulators retained per target gene; also sets maxModSize for R²
xaxisStepSize    = Int(meanEdgesPerGene / 4)   # x-axis tick spacing on model_size.pdf
correlationWeight = 1          # weight of partial correlation in edge ranking; 0 = stability score only
lambdaBias       = [0.5]       # prior penalty weight(s): 0 = ignore prior, 1 = full prior penalty

# ── Steps 1–3: load data once ─────────────────────────────────────────────────
println("Loading expression data and prior...")
data                 = loadData(normGeneExprFile, targGeneFile, potRegFile)
priorData, mergedTFs = loadPrior(data, priorFile)
estimateTFA(priorData, data)

# ── Main loop: leave-out sets × predictor modes ───────────────────────────────
for lind in 1:size(loInfo, 1)
    leaveOutFile  = loInfo[lind][1]
    leaveOutLabel = loInfo[lind][2]
    println("\n=== Leave-out: $leaveOutLabel ===")

    for (tfaMode, modeLabel) in zip(tfaOptions, modeLabels)
        println("  Mode: $modeLabel")

        dirOut = joinpath(outputDir, leaveOutLabel, modeLabel)
        mkpath(dirOut)

        # ── Step 4: build GRN with leave-out samples excluded ─────────────────
        # buildNetwork saves instabOutMat.jld (serialised GrnData) to dirOut
        # when outputDir is provided. trainInds stored inside GrnData reflects
        # the leave-out split and is used by calcR2predFromStabilities.
        buildGrn = buildNetwork(data, priorData;
                                tfaMode              = tfaMode,
                                priorFilePenalties   = priorFilePenalties,
                                lambdaBias           = lambdaBias,
                                totSS                = totSS,
                                bstarsTotSS          = bstarsTotSS,
                                subsampleFrac        = subsampleFrac,
                                minLambda            = lambdaMin,
                                maxLambda            = lambdaMax,
                                totLambdasBstars     = totLambdasBstars,
                                totLambdas           = totLambdas,
                                targetInstability    = targetInstability,
                                meanEdgesPerGene     = meanEdgesPerGene,
                                correlationWeight    = correlationWeight,
                                instabilityLevel     = instabilityLevel,
                                modelSelection       = modelSelection,
                                leaveOutSampleList   = leaveOutFile,
                                outputDir            = dirOut)

        # ── Step 5: reload GrnData written by bstartsEstimateInstability ──────
        # buildNetwork saves instabOutMat.jld; load it to get trainInds,
        # responseMat, predictorMat, subsamps needed for R² computation.
        grnData = load_object(joinpath(dirOut, "instabOutMat.jld"))

        # ── Step 6: compute R²_pred across model sizes ────────────────────────
        println("  Computing R² prediction performance...")
        calcR2predFromStabilities(grnData, buildGrn, data,
                                  meanEdgesPerGene, dirOut;
                                  xaxisStepSize = xaxisStepSize)
        println("  Outputs → $dirOut")
    end
end

println("\nDone. All R² outputs saved under: $outputDir")
