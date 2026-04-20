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

normGeneExprFile = "/path/to/expression.arrow"     # genes × samples (Arrow or TSV)
targGeneFile     = "/path/to/target_genes.txt"
potRegFile       = "/path/to/pot_regs.txt"
priorFile        = "/path/to/prior.tsv"
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

totSS            = 80
bstarsTotSS      = 5
subsampleFrac    = 0.63
lambdaBias       = [0.5]
targetInstability = 0.05
lambdaMin        = 0.01
lambdaMax        = 1.0
totLambdasBstars = 20
totLambdas       = 40
meanEdgesPerGene = 20
xaxisStepSize    = Int(meanEdgesPerGene / 4)
correlationWeight = 1
instabilityLevel = "Network"

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
