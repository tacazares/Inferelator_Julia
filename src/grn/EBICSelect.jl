# =============================================================================
#  EBICSelect.jl — EBIC and bEBIC model selection for LASSO-based GRN inference
#
#  Provides two alternative lambda selection methods to bStARS:
#
#  EBIC  (Extended BIC): fits LASSO on the full dataset once per target gene
#        and selects the lambda that minimises the Extended Bayesian Information
#        Criterion. Fast (no subsampling), but less robust to noise.
#
#  bEBIC (bootstrap EBIC): fits LASSO on each subsample independently, selects
#        the EBIC-optimal lambda per subsample, then uses the median lambda
#        across subsamples for a final full-data fit. Selection frequency across
#        subsamples at each subsample's EBIC-optimal lambda is used as the edge
#        confidence score — analogous to the bStARS stability score, and
#        compatible with rankEdges!.
#
#  Both methods use the same prior-weighted adaptive LASSO penalties
#  (penaltyMat) as bStARS. The output (buildGrn.networkStability, buildGrn.betas,
#  buildGrn.lambda) is compatible with rankEdges! so the rest of the pipeline
#  is unchanged.
#
#  Reference:
#    Chen & Chen (2008) "Extended Bayesian Information Criteria for Model
#    Selection with Large Model Spaces." Biometrika 95(3):759-771.
# =============================================================================


"""
    computeEBIC(betas, responses, predictors, n, p; gamma=0.5) -> Vector{Float64}

Compute the Extended BIC score at each point along a LASSO solution path.

# Arguments
- `betas`      : Coefficient matrix (n_predictors x n_lambdas) from GLMNet
- `responses`  : Response vector (n samples)
- `predictors` : Predictor matrix (n samples x p_predictors)
- `n`          : Number of samples used for fitting
- `p`          : Number of candidate predictors (before penalty-induced sparsity)
- `gamma`      : EBIC hyperparameter controlling the sparsity penalty (default 0.5).
                 gamma=0 reduces to standard BIC.
                 gamma=0.5 is recommended for high-dimensional sparse settings (GRN).
                 gamma=1 maximises the sparsity penalty.

# Returns
`Vector{Float64}` of EBIC scores, one per lambda on the solution path. Lower = better.

# Details
EBIC(lambda) = n * log(RSS/n) + k * log(n) + 2 * gamma * log C(p, k)

where:
- RSS = residual sum of squares at lambda
- k   = number of non-zero coefficients at lambda
- C(p,k) = binomial coefficient, computed via lgamma for numerical stability.

The first term measures fit quality; the second penalises model complexity
(standard BIC); the third adds a graph-selection penalty that grows with the
number of possible models of size k from p candidates.
"""
function computeEBIC(
    betas::AbstractMatrix,
    responses::AbstractVector,
    predictors::AbstractMatrix,
    n::Int, p::Int;
    gamma::Float64 = 0.5
)::Vector{Float64}

    nLambdas   = size(betas, 2)
    ebicScores = Vector{Float64}(undef, nLambdas)

    for li in 1:nLambdas
        coefs = vec(betas[:, li])
        k     = sum(coefs .!= 0)

        # Residual sum of squares
        resid = responses .- predictors * coefs
        rss   = sum(resid .^ 2)
        rss   = max(rss, 1e-10)   # guard against exact-zero RSS (numerical edge case)

        # Gaussian log-likelihood term: -2 * logLik = n * log(RSS / n)
        logLikTerm = n * log(rss / n)

        # log C(p, k) via lgamma for numerical stability
        logBinom = (k == 0 || k >= p) ? 0.0 :
                   lgamma(p + 1) - lgamma(k + 1) - lgamma(p - k + 1)

        ebicScores[li] = logLikTerm + k * log(n) + 2 * gamma * logBinom
    end

    return ebicScores
end


"""
    ebicSelect!(grnData, buildGrn; kwargs...)

Select the regularisation lambda per target gene using a single full-data LASSO
fit evaluated by the Extended BIC (EBIC). This is the fastest model selection
alternative to bStARS — no subsampling is performed.

Edge confidence scores in `buildGrn.networkStability` are set to the absolute
coefficient values at the EBIC-optimal lambda, making the output directly
compatible with `rankEdges!`.

# Arguments
- `grnData`    : `GrnData` with `penaltyMat`, `predictorMat`, and `responseMat` populated
- `buildGrn`   : `BuildGrn` struct — populated in-place with `networkStability`,
                 `betas`, and `lambda` fields
- `gamma`      : EBIC hyperparameter (default 0.5; see `computeEBIC` for details)
- `minLambda`  : Lower bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `maxLambda`  : Upper bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `totLambdas` : Number of grid points when `minLambda`/`maxLambda` are user-supplied
                 (default 100). Ignored when either bound is `nothing`.
- `zScoreLASSO`: Z-score predictors and (optionally) response before fitting (default true)

# Notes
- Supply `minLambda`/`maxLambda` if you want to constrain the search range
  (e.g., to a region informed by a prior bStARS run). Leave as `nothing` for
  an unconstrained, data-driven solution path.
- Unlike bStARS, EBIC produces no instability-based confidence: the edge
  confidence is the absolute LASSO coefficient, which encodes effect size
  rather than reproducibility. Use bEBIC (`bebicSelect!`) if you prefer a
  selection-frequency-based confidence score.
"""
function ebicSelect!(
    grnData::GrnData,
    buildGrn::BuildGrn;
    gamma::Float64                     = 0.5,
    minLambda::Union{Float64, Nothing} = nothing,
    maxLambda::Union{Float64, Nothing} = nothing,
    totLambdas::Int                    = 100,
    zScoreLASSO::Bool                  = true
)
    totResponses  = size(grnData.responseMat, 1)
    totPredictors = size(grnData.predictorMat, 1)
    nSamples      = size(grnData.responseMat, 2)

    # Pre-compute finite-predictor indices (Inf penalty = illegal interaction)
    responsePredInds = Vector{Vector{Int}}(undef, totResponses)
    for res in 1:totResponses
        responsePredInds[res] = findall(x -> x != Inf, grnData.penaltyMat[res, :])
    end

    networkStability = zeros(totResponses, totPredictors)
    betas_out        = zeros(totResponses, totPredictors)
    ebicLambdas      = zeros(totResponses)

    @info "EBIC lambda selection — fitting full-data LASSO per gene"
    Threads.@threads for res in ProgressBar(1:totResponses)
        predInds      = responsePredInds[res]
        currPredNum   = length(predInds)
        penaltyFactor = grnData.penaltyMat[res, predInds]
        p             = currPredNum

        # Z-score predictors (all samples)
        dt        = fit(ZScoreTransform, grnData.predictorMat[predInds, :], dims=2)
        currPreds = transpose(StatsBase.transform(dt, grnData.predictorMat[predInds, :]))

        if zScoreLASSO
            dtR           = fit(ZScoreTransform, grnData.responseMat[res, :], dims=1)
            currResponses = StatsBase.transform(dtR, grnData.responseMat[res, :])
        else
            currResponses = grnData.responseMat[res, :]
        end

        # Fit LASSO: use user-supplied range or let GLMNet choose
        if minLambda !== nothing && maxLambda !== nothing
            lambdaGrid = reverse(collect(range(minLambda, maxLambda, length = totLambdas)))
            lsoln = glmnet(currPreds, vec(currResponses);
                           penalty_factor = penaltyFactor,
                           lambda         = lambdaGrid,
                           alpha          = 1.0)
        else
            lsoln = glmnet(currPreds, vec(currResponses);
                           penalty_factor = penaltyFactor,
                           alpha          = 1.0)
        end

        # Select lambda at minimum EBIC
        ebicScores = computeEBIC(lsoln.betas, vec(currResponses), Matrix(currPreds),
                                 nSamples, p; gamma = gamma)
        bestInd    = argmin(ebicScores)
        bestCoefs  = vec(lsoln.betas[:, bestInd])

        ebicLambdas[res]                    = lsoln.lambda[bestInd]
        networkStability[res, predInds]     = abs.(bestCoefs)   # confidence = |coef|
        betas_out[res, predInds]            = bestCoefs
    end

    # Zero-out Inf entries (illegal TF-gene pairs)
    nsVec             = networkStability[:]
    nsVec[isinf.(nsVec)] .= 0.0
    networkStability  = reshape(nsVec, totResponses, totPredictors)

    buildGrn.networkStability = networkStability
    buildGrn.betas            = betas_out
    buildGrn.lambda           = ebicLambdas
    grnData.ebicLambdas       = ebicLambdas

    nSelected = sum(ebicLambdas .> 0)
    @info "EBIC complete" genesProcessed=totResponses medianLambda=median(ebicLambdas)
end


"""
    bebicSelect!(grnData, buildGrn; kwargs...)

Select the regularisation lambda per target gene using bootstrap EBIC (bEBIC).

For each target gene, LASSO is fitted on each subsample and the EBIC-optimal
lambda is recorded. The median lambda across subsamples becomes the gene's final
lambda. Selection frequency across subsamples at each subsample's EBIC-optimal
lambda is used as the edge confidence score — analogous to the bStARS selection
frequency theta, and compatible with `rankEdges!`.

A final LASSO fit on the full dataset at the median lambda extracts the signed
coefficients stored in `buildGrn.betas`.

# Arguments
- `grnData`    : `GrnData` with `penaltyMat`, `predictorMat`, `responseMat`, and
                 `subsamps` populated (call `constructSubsamples` first)
- `buildGrn`   : `BuildGrn` struct — populated in-place
- `gamma`      : EBIC hyperparameter (default 0.5; see `computeEBIC` for details)
- `minLambda`  : Lower bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `maxLambda`  : Upper bound of the lambda grid. `nothing` lets GLMNet choose.
- `totLambdas` : Number of grid points when `minLambda`/`maxLambda` are supplied
                 (default 100).
- `zScoreLASSO`: Z-score predictors before each fit (default true)

# Notes
- Unlike pure EBIC, bEBIC produces a selection-frequency confidence score
  (0 to 1), matching the output semantics of bStARS. This makes bEBIC the
  preferred choice when you want results directly comparable to bStARS.
- `constructSubsamples` must be called before `bebicSelect!` to populate
  `grnData.subsamps`.
- The subsamples used are the same as those passed to bStARS, controlled by
  `totSS` and `subsampleFrac` in `constructSubsamples`.
"""
function bebicSelect!(
    grnData::GrnData,
    buildGrn::BuildGrn;
    gamma::Float64                     = 0.5,
    minLambda::Union{Float64, Nothing} = nothing,
    maxLambda::Union{Float64, Nothing} = nothing,
    totLambdas::Int                    = 100,
    zScoreLASSO::Bool                  = true
)
    totResponses  = size(grnData.responseMat, 1)
    totPredictors = size(grnData.predictorMat, 1)
    nSamples      = size(grnData.responseMat, 2)
    totSS         = size(grnData.subsamps, 1)

    if totSS == 0
        error("bEBIC requires subsamples. Call constructSubsamples before bebicSelect!.")
    end

    # Pre-compute finite-predictor indices
    responsePredInds = Vector{Vector{Int}}(undef, totResponses)
    for res in 1:totResponses
        responsePredInds[res] = findall(x -> x != Inf, grnData.penaltyMat[res, :])
    end

    networkStability = zeros(totResponses, totPredictors)
    betas_out        = zeros(totResponses, totPredictors)
    ebicLambdas      = zeros(totResponses)

    @info "bEBIC lambda selection — fitting LASSO on $totSS subsamples per gene"
    Threads.@threads for res in ProgressBar(1:totResponses)
        predInds      = responsePredInds[res]
        currPredNum   = length(predInds)
        penaltyFactor = grnData.penaltyMat[res, predInds]
        p             = currPredNum

        bestLambdasSS  = zeros(totSS)
        ssSelections   = zeros(totSS, currPredNum)   # binary selection at each subsample's EBIC-optimal lambda

        # -- Subsample loop: find EBIC-optimal lambda per subsample --
        for ss in 1:totSS
            subsamp = grnData.subsamps[ss, :]
            nSS     = length(subsamp)

            dt        = fit(ZScoreTransform, grnData.predictorMat[predInds, subsamp], dims=2)
            currPreds = transpose(StatsBase.transform(dt, grnData.predictorMat[predInds, subsamp]))

            if zScoreLASSO
                dtR           = fit(ZScoreTransform, grnData.responseMat[res, subsamp], dims=1)
                currResponses = StatsBase.transform(dtR, grnData.responseMat[res, subsamp])
            else
                currResponses = grnData.responseMat[res, subsamp]
            end

            if minLambda !== nothing && maxLambda !== nothing
                lambdaGrid = reverse(collect(range(minLambda, maxLambda, length = totLambdas)))
                lsoln = glmnet(currPreds, vec(currResponses);
                               penalty_factor = penaltyFactor,
                               lambda         = lambdaGrid,
                               alpha          = 1.0)
            else
                lsoln = glmnet(currPreds, vec(currResponses);
                               penalty_factor = penaltyFactor,
                               alpha          = 1.0)
            end

            ebicScores         = computeEBIC(lsoln.betas, vec(currResponses),
                                             Matrix(currPreds), nSS, p; gamma = gamma)
            bestInd            = argmin(ebicScores)
            bestLambdasSS[ss]  = lsoln.lambda[bestInd]
            # Binary selection at this subsample's EBIC-optimal lambda
            ssSelections[ss, :] = Float64.(vec(lsoln.betas[:, bestInd]) .!= 0)
        end

        # Median lambda across subsamples
        finalLambda = median(bestLambdasSS)

        # Selection frequency = average binary selection across subsamples
        # (analogous to theta in bStARS — used as edge confidence score)
        selectionFreq = vec(mean(ssSelections, dims=1))

        # -- Final fit on full data at median lambda --
        dtFull    = fit(ZScoreTransform, grnData.predictorMat[predInds, :], dims=2)
        fullPreds = transpose(StatsBase.transform(dtFull, grnData.predictorMat[predInds, :]))

        if zScoreLASSO
            dtRFull       = fit(ZScoreTransform, grnData.responseMat[res, :], dims=1)
            fullResponses = StatsBase.transform(dtRFull, grnData.responseMat[res, :])
        else
            fullResponses = grnData.responseMat[res, :]
        end

        lsolnFull = glmnet(fullPreds, vec(fullResponses);
                           penalty_factor = penaltyFactor,
                           lambda         = [finalLambda],
                           alpha          = 1.0)
        bestCoefs = vec(lsolnFull.betas[:, 1])

        ebicLambdas[res]                = finalLambda
        networkStability[res, predInds] = selectionFreq   # confidence = selection frequency
        betas_out[res, predInds]        = bestCoefs
    end

    # Zero-out Inf entries
    nsVec                = networkStability[:]
    nsVec[isinf.(nsVec)] .= 0.0
    networkStability     = reshape(nsVec, totResponses, totPredictors)

    buildGrn.networkStability = networkStability
    buildGrn.betas            = betas_out
    buildGrn.lambda           = ebicLambdas
    grnData.ebicLambdas       = ebicLambdas

    @info "bEBIC complete" genesProcessed=totResponses medianLambda=median(ebicLambdas) totSubsamples=totSS
end
