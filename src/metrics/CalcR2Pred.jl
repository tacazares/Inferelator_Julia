"""
    calcR2predFromStabilities(grnData, buildGrn, expressionData, maxModSize, outputDir; xaxisStepSize=5)

Estimate out-of-sample predictive performance (R²_pred) as a function of network
model size, using confidence-weighted edge selection on held-out samples.

Compatible with bStARS and bEBIC model selection. For each model size m
(1 → maxModSize TFs per gene) a confidence cutoff is derived from
buildGrn.networkStability so that on average m edges per gene are retained.
An OLS model is fit on the training samples (grnData.trainInds) and evaluated on
the held-out test samples (all samples NOT in trainInds). Four metrics are reported
at each model size:

- R²_pred   : out-of-sample R², all target gene models (including null/zero models)
- R²_fit    : in-sample R², all target gene models
- R²_predNz : R²_pred restricted to non-zero models (at least one edge selected)
- R²_fitNz  : R²_fit restricted to non-zero models

Requires a leave-out sample set; when trainInds == 1:totSamples (no leave-out),
the test set is empty and R²_pred cannot be computed — an error is thrown.

# Arguments
- `grnData`        : GrnData after constructSubsamples + instability estimation.
                     Uses trainInds, responseMat, predictorMat, subsamps.
                     Compatible with bStARS and bEBIC (both call constructSubsamples).
- `buildGrn`       : BuildGrn after model selection + rankEdges!.
                     Uses networkStability:
                       bStARS → raw subsample counts (0 to totSS)
                       bEBIC  → selection frequencies (0 to 1) — auto-detected
- `expressionData` : GeneExpressionData, provides cellLabels for total sample count.
- `maxModSize`     : Maximum TFs per gene to sweep; 1:maxModSize model sizes tested.
- `outputDir`      : Directory for r2_pred.pdf, model_size.pdf, and r2pred.jld2.
- `xaxisStepSize`  : X-axis tick spacing for model_size.pdf (default 5).

# Outputs (written to outputDir)
- `r2_pred.pdf`    : R²_pred and R²_fit vs selection frequency cutoff.
- `model_size.pdf` : R²_pred vs average number of TFs per gene.
- `r2pred.jld2`    : Saved results dict (instabRangeNet, numParams, r2_pred, …).
"""
function calcR2predFromStabilities(
    grnData::GrnData,
    buildGrn::BuildGrn,
    expressionData::GeneExpressionData,
    maxModSize::Int,
    outputDir::String;
    xaxisStepSize::Int = 5
)
    targGenes     = grnData.targGenes
    allPredictors = grnData.allPredictors
    conditionsc   = expressionData.cellLabels
    allStabsTest  = buildGrn.networkStability      # targets × predictors: counts (bStARS) or frequencies (bEBIC)
    trainInds     = grnData.trainInds
    responseMat   = grnData.responseMat
    predictorMat  = grnData.predictorMat
    totSS         = size(grnData.subsamps, 1)       # derived from stored subsample matrix

    # Auto-detect confidence scale:
    # bStARS stores raw subsample counts (0 to totSS, max > 1)
    # bEBIC  stores selection frequencies (0 to 1, max <= 1)
    stabNormFactor = (totSS > 0 && maximum(allStabsTest) > 1.0) ? Float64(totSS) : 1.0

    totPredTargs = length(targGenes)
    totConds     = length(conditionsc)
    totPreds     = length(allPredictors)
    modSizes     = collect(1:maxModSize)
    totModSizes  = length(modSizes)

    testInds = setdiff(1:totConds, trainInds)
    if isempty(testInds)
        error("calcR2predFromStabilities: testInds is empty. " *
              "A leave-out sample set is required for R²_pred computation. " *
              "Pass leaveOutSampleList to buildNetwork (or constructSubsamples).")
    end

    # ── Stability cutoffs corresponding to each model size ─────────────────────
    # Sort edges by descending stability count and find the count at which
    # floor(totPredTargs * m) edges are retained — this is the cutoff for model size m.
    instabRangeNet   = zeros(totModSizes)
    inds             = reverse(sortperm(vec(allStabsTest)))
    stabOrd          = vec(allStabsTest)[inds]
    stabOrd          = stabOrd[findall(!iszero, stabOrd)]
    totInts          = length(stabOrd)
    subsampleCutoffs = Float64.(modSizes)

    for ms = 1:totModSizes
        currNetSize = floor(Int, totPredTargs * modSizes[ms])
        if currNetSize <= totInts
            subsampleCutoffs[ms] = stabOrd[currNetSize]
            currTheta = subsampleCutoffs[ms] / stabNormFactor
            instabRangeNet[ms] = 2 * currTheta * (1 - currTheta)
        end
        # if currNetSize > totInts the cutoff stays at its initialised value (Float64(ms));
        # the corresponding model size is larger than the number of non-zero edges observed
    end

    # ── Allocate SSE accumulators ───────────────────────────────────────────────
    kfoldCv = 1   # single leave-out split
    totTrainConds = length(trainInds)
    totTestConds  = length(testInds)

    SSE_predMat          = zeros(totModSizes, totConds,    totPredTargs)
    SSE_fitMat           = zeros(totModSizes, totConds,    totPredTargs)
    SSE_fitMatMean       = zeros(totConds,    totPredTargs)
    SSE_predMatMean      = zeros(totConds,    totPredTargs)
    SSE_predMatNz        = zeros(totModSizes)
    SSE_predMatMeanNz    = zeros(totModSizes)
    SSE_fitMatNz         = zeros(totModSizes)
    SSE_fitMatMeanNz     = zeros(totModSizes)
    r2_pred_byTarg       = zeros(totModSizes, totPredTargs)
    r2_fit_byTargs       = zeros(totModSizes, totPredTargs)
    r2_pred_byCond       = zeros(totModSizes, totConds)
    r2_fit_byCond        = zeros(totModSizes, totConds)
    r2_pred              = zeros(totModSizes)
    r2_fit               = zeros(totModSizes)
    r2_predNz            = zeros(totModSizes)
    r2_fitNz             = zeros(totModSizes)
    numParamsMat         = zeros(totModSizes, kfoldCv, totPredTargs)
    numParams            = zeros(totModSizes, totPredTargs)

    # ── Single leave-out fold ───────────────────────────────────────────────────
    targGenesTrain = responseMat[:, trainInds]
    targGenesTest  = responseMat[:, testInds]
    tfaTrain       = predictorMat[:, trainInds]
    tfaTest        = predictorMat[:, testInds]

    # Mean-centre targets using training mean; z-score TFA using training stats
    meanTargGene    = mean(targGenesTrain, dims=2)
    cTargGenesTrain = targGenesTrain .- meanTargGene
    cTargGenesTest  = targGenesTest  .- meanTargGene
    meanTfa = mean(tfaTrain, dims=2)
    stdTfa  = std(tfaTrain,  dims=2)
    zTfaTrain = (tfaTrain .- meanTfa) ./ stdTfa
    zTfaTest  = (tfaTest  .- meanTfa) ./ stdTfa

    # Baseline SSE (null model = training mean)
    SSE_predMatMean[testInds,  :] = cTargGenesTest'.^2
    SSE_fitMatMean[trainInds,  :] = cTargGenesTrain'.^2

    for targ = 1:totPredTargs
        currTargValsTrain = cTargGenesTrain[targ, :]
        currTargValsTest  = cTargGenesTest[targ,  :]
        currSubsamples    = allStabsTest[targ, :]

        for iind = 1:totModSizes
            currCut   = subsampleCutoffs[iind]
            intInds   = findall(x -> x >= currCut, currSubsamples)
            totParams = length(intInds)
            numParamsMat[iind, 1, targ] = totParams
            numParams[iind, targ]      += totParams

            if totParams > 0
                currPredValsTrain = zTfaTrain[intInds, :]'   # totTrainConds × totParams
                currPredValsTest  = zTfaTest[intInds,  :]'   # totTestConds  × totParams
                coefs      = currPredValsTrain \ currTargValsTrain
                trainPreds = currPredValsTrain * coefs        # totTrainConds vector
                testPreds  = currPredValsTest  * coefs        # totTestConds  vector

                SSE_predMatNz[iind]     += sum((testPreds  .- currTargValsTest).^2)
                SSE_predMatMeanNz[iind] += sum(SSE_predMatMean[testInds,  targ])
                SSE_fitMatNz[iind]      += sum((trainPreds .- currTargValsTrain).^2)
                SSE_fitMatMeanNz[iind]  += sum(SSE_fitMatMean[trainInds, targ])
            else
                trainPreds = zeros(totTrainConds)
                testPreds  = zeros(totTestConds)
            end

            SSE_predMat[iind, testInds,  targ]  = (testPreds  .- currTargValsTest).^2
            SSE_fitMat[iind,  trainInds, targ] .+= (trainPreds .- currTargValsTrain).^2
        end
    end

    # ── R² ─────────────────────────────────────────────────────────────────────
    r2_predNz = 1 .- SSE_predMatNz ./ SSE_predMatMeanNz
    r2_fitNz  = 1 .- SSE_fitMatNz  ./ SSE_fitMatMeanNz

    for bindd = 1:totModSizes
        currSSE_preds = SSE_predMat[bindd, :, :]   # totConds × totPredTargs
        currSSE_fit   = SSE_fitMat[bindd,  :, :]
        r2_pred[bindd] = 1 - sum(currSSE_preds) / sum(SSE_predMatMean)
        r2_fit[bindd]  = 1 - sum(currSSE_fit)   / sum(SSE_fitMatMean)
        r2_pred_byTarg[bindd, :] = vec(1 .- sum(currSSE_preds, dims=1) ./ sum(SSE_predMatMean, dims=1))
        r2_fit_byTargs[bindd, :] = vec(1 .- sum(currSSE_fit,  dims=1) ./ sum(SSE_fitMatMean,  dims=1))
        r2_pred_byCond[bindd, :] = vec(1 .- sum(currSSE_preds', dims=1) ./ sum(SSE_predMatMean', dims=1))
        r2_fit_byCond[bindd,  :] = vec(1 .- sum(currSSE_fit',  dims=1) ./ sum(SSE_fitMatMean',  dims=1))
    end

    # ── Plots ───────────────────────────────────────────────────────────────────
    mkpath(outputDir)
    axis_title_size = 15
    tick_label_size = 13

    # Plot 1 — R² vs instability cutoff
    figure(figsize=(3.6, 4))
    plot(instabRangeNet, r2_pred,   "b-",  label="Predicted")
    plot(instabRangeNet, r2_fit,    "y-",  label="Fit")
    plot(instabRangeNet, r2_predNz, "g--", label="nzPredicted")
    plot(instabRangeNet, r2_fitNz,  "r--", label="nzFit")
    xlabel("Selection frequency cutoff", fontsize=axis_title_size)
    ylabel(L"R^2\ \left(1 - \frac{SSE_{model}}{SSE_{mean}}\right)", fontsize=axis_title_size)
    legend()
    grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")
    grid(true, which="minor", linestyle=":",  linewidth=0.5,  color="lightgray")
    minorticks_on()
    tick_params(axis="both", which="major", labelsize=tick_label_size)
    ylim(0, 1)
    gca().set_yticks(0:0.2:1)
    tight_layout()
    savefig(joinpath(outputDir, "r2_pred.pdf"), dpi=600)
    close()

    # Plot 2 — R² vs average TFs per gene
    totEdges = vec(sum(numParams, dims=2) ./ totPredTargs)
    figure(figsize=(3.6, 4))
    plot(totEdges,                         r2_pred,                             "k-")
    plot(vec(median(numParams, dims=2)),   vec(median(r2_pred_byTarg, dims=2)), "b-")
    plot(vec(mean(numParams,   dims=2)),   vec(mean(r2_pred_byTarg,   dims=2)), "r-")
    xlabel("# of Predictors per Gene", fontsize=axis_title_size)
    ylabel(L"R^2\ \left(1 - \frac{SSE_{model}}{SSE_{mean}}\right)", fontsize=axis_title_size)
    legend([L"Overall\ R^2_{pred}", L"Gene\ Median\ R^2_{pred}", L"Gene\ Average\ R^2_{pred}"])
    grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")
    grid(true, which="minor", linestyle=":",  linewidth=0.5,  color="lightgray")
    minorticks_on()
    tick_params(axis="both", which="major", labelsize=tick_label_size)
    ylim(0, 1)
    gca().set_yticks(0:0.2:1)
    gca().set_xticks(0:xaxisStepSize:maxModSize)
    tight_layout()
    savefig(joinpath(outputDir, "model_size.pdf"), dpi=600)
    close()

    # ── Save ────────────────────────────────────────────────────────────────────
    r2OutMat = joinpath(outputDir, "r2pred.jld2")
    jldsave(r2OutMat;
            instabRangeNet, numParams, r2_pred, r2_fit, r2_predNz, r2_fitNz,
            totEdges, r2_pred_byTarg, r2_pred_byCond)
    @info "R² results saved → $r2OutMat"
end
