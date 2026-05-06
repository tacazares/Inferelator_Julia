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
LASSO coefficient at the chosen lambda (effect-size scale, not reproducibility).

# Saved diagnostics
After the run, `grnData` contains:
- `ebicLambdas`  : chosen lambda per gene (also mirrored in `buildGrn.lambda`)
- `ebicNonzero`  : number of non-zero coefficients at the chosen lambda per gene

Use `writeEBICLambdaTable!` to write these to a TSV file.

# Arguments
- `grnData`    : `GrnData` with `penaltyMat`, `predictorMat`, and `responseMat` populated
- `buildGrn`   : `BuildGrn` struct — populated in-place with `networkStability`,
                 `betas`, and `lambda` fields
- `gamma`      : EBIC hyperparameter (default 0.5; see `computeEBIC` for details)
- `minLambda`  : Lower bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `maxLambda`  : Upper bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `totLambdas` : Number of grid points when `minLambda`/`maxLambda` are
                 user-supplied (default 100). Ignored when either bound is
                 `nothing`. Grid is linear by default; see inline comment in
                 source for a log-spaced alternative.
- `zScoreLASSO`: Z-score predictors and (optionally) response before fitting (default true)

# Notes
- Supply `minLambda`/`maxLambda` if you want to constrain the search range
  (e.g., to a region informed by a prior bStARS run). Leave as `nothing` for
  an unconstrained, data-driven solution path.
- Unlike bStARS, EBIC produces no instability-based confidence: the edge
  confidence is the absolute LASSO coefficient, which encodes effect size
  rather than reproducibility. Use bEBIC (`bebicEstimateLambdas!` +
  `bebicFinalFit!`) if you prefer a selection-frequency-based confidence score.
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

    # Warn when response is not z-scored: lambda_max scales with response magnitude,
    # so GLMNet's auto-chosen grid will be much larger than the bStARS range.
    # EBIC may then select an overly sparse (or null) model.
    # Fix: set zScoreLASSO=true, or supply minLambda/maxLambda explicitly.
    if !zScoreLASSO && isnothing(minLambda)
        @warn "ebicSelect!: zScoreLASSO=false and no minLambda/maxLambda supplied. " *
              "The response is on its raw scale while predictors are z-scored, so " *
              "lambda_max ∝ response magnitude and GLMNet's auto grid may be very large. " *
              "EBIC is likely to select an overly sparse model. " *
              "Recommended fix: set zScoreLASSO=true, or pass minLambda/maxLambda."
    end

    networkStability  = zeros(totResponses, totPredictors)
    betas_out         = zeros(totResponses, totPredictors)
    ebicLambdas       = zeros(totResponses)
    nNonzeroVec       = zeros(Int, totResponses)   # non-zero coefficients at chosen λ per gene
    # Boundary hit trackers — each thread writes to its own res index, no race condition.
    # upperBoundaryHits[res]=1 → bestInd==1 (most regularised end) → maxLambda too small.
    # lowerBoundaryHits[res]=1 → bestInd==end (least regularised end) → minLambda too large.
    upperBoundaryHits = zeros(Int, totResponses)
    lowerBoundaryHits = zeros(Int, totResponses)

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
            # Log-spaced alternative — concentrates steps in the low-lambda region
            # where EBIC typically finds its minimum, matching GLMNet's auto path.
            # Swap in if linear grid undersamples the optimal region:
            # lambdaGrid = reverse(exp.(range(log(minLambda), log(maxLambda), length = totLambdas)))
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
        bestInd   = argmin(ebicScores)
        bestCoefs = vec(lsoln.betas[:, bestInd])

        # Boundary detection (only meaningful on a user-supplied grid)
        if minLambda !== nothing && maxLambda !== nothing
            if bestInd == 1
                upperBoundaryHits[res] = 1   # optimum at most-regularised end → maxLambda too small
            elseif bestInd == length(lsoln.lambda)
                lowerBoundaryHits[res] = 1   # optimum at least-regularised end → minLambda too large
            end
        end

        ebicLambdas[res]                = lsoln.lambda[bestInd]
        nNonzeroVec[res]                = sum(bestCoefs .!= 0)   # model size at chosen λ
        networkStability[res, predInds] = abs.(bestCoefs)        # confidence = |coef|
        betas_out[res, predInds]        = bestCoefs
    end

    # Boundary warnings — fired once after the loop with aggregate counts
    if minLambda !== nothing && maxLambda !== nothing
        nUpper = sum(upperBoundaryHits)
        nLower = sum(lowerBoundaryHits)
        if nUpper > 0
            @warn "ebicSelect!: EBIC minimum at upper boundary (maxLambda=$maxLambda) " *
                  "for $nUpper / $totResponses genes. True optimum may lie above maxLambda — " *
                  "consider increasing maxLambda or setting it to nothing for automatic path selection."
        end
        if nLower > 0
            @warn "ebicSelect!: EBIC minimum at lower boundary (minLambda=$minLambda) " *
                  "for $nLower / $totResponses genes. True optimum may lie below minLambda — " *
                  "consider decreasing minLambda or setting it to nothing for automatic path selection."
        end
    end

    # Zero-out Inf entries (illegal TF-gene pairs)
    nsVec             = networkStability[:]
    nsVec[isinf.(nsVec)] .= 0.0
    networkStability  = reshape(nsVec, totResponses, totPredictors)

    buildGrn.networkStability = networkStability
    buildGrn.betas            = betas_out
    buildGrn.lambda           = ebicLambdas
    grnData.ebicLambdas       = ebicLambdas
    grnData.ebicNonzero       = nNonzeroVec

    @info "EBIC complete" genesProcessed=totResponses medianLambda=median(ebicLambdas) medianNonzero=median(nNonzeroVec)
end


"""
    bebicEstimateLambdas!(grnData, buildGrn; kwargs...)

For each target gene, fit LASSO on every subsample and select the EBIC-optimal
lambda per subsample. Records the per-subsample lambda distribution in
`grnData.lambdaSS` and the selection counts in `buildGrn.networkStability`.

This is the expensive, parallelised stage of the bEBIC pipeline. Call
`bebicFinalFit!` afterwards to aggregate the per-subsample lambdas and produce
the final full-data fit.

# Arguments
- `grnData`    : `GrnData` with `penaltyMat`, `predictorMat`, `responseMat`, and
                 `subsamps` populated (call `constructSubsamples` first)
- `buildGrn`   : `BuildGrn` struct — `networkStability` populated in-place
- `gamma`      : EBIC hyperparameter (default 0.5; see `computeEBIC`)
- `minLambda`  : Lower bound of the lambda grid. `nothing` (default) lets GLMNet
                 choose its own solution path automatically (recommended).
- `maxLambda`  : Upper bound of the lambda grid. `nothing` lets GLMNet choose.
- `totLambdas` : Number of grid points when `minLambda`/`maxLambda` are supplied
                 (default 100). Ignored when either bound is `nothing`. Grid is
                 linear by default; see inline comment in source for a log-spaced
                 alternative.
- `zScoreLASSO`: Z-score predictors before each subsample fit (default true)

# Fills
- `grnData.lambdaSS`          : `totResponses × totSS` matrix — the EBIC-optimal
                                 lambda selected on each subsample for each gene.
                                 Row `i` is the per-subsample lambda distribution
                                 for gene `i`; high std indicates unstable selection.
- `buildGrn.networkStability` : selection counts (0–totSS) per gene × predictor pair,
                                 analogous to bStARS's stability score

# Notes
- `constructSubsamples` must be called first to populate `grnData.subsamps`.
- `gamma` is applied during EBIC scoring here; changing it requires re-running
  this function (not just `bebicFinalFit!`).
- Unlike pure EBIC, bEBIC produces a selection-count confidence score (0–totSS)
  matching the output scale of bStARS.
"""
function bebicEstimateLambdas!(
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
    totSS         = size(grnData.subsamps, 1)

    if totSS == 0
        error("bEBIC requires subsamples. Call constructSubsamples before bebicEstimateLambdas!.")
    end

    # Pre-compute finite-predictor indices
    responsePredInds = Vector{Vector{Int}}(undef, totResponses)
    for res in 1:totResponses
        responsePredInds[res] = findall(x -> x != Inf, grnData.penaltyMat[res, :])
    end

    networkStability  = zeros(totResponses, totPredictors)
    # Per-subsample EBIC-optimal lambda matrix: row = gene, col = subsample.
    # Each cell is the lambda that minimised EBIC on that subsample for that gene.
    # std(row) measures how stable the EBIC-optimal lambda is across subsamples.
    lambdaSSMat       = zeros(totResponses, totSS)
    # Boundary hit counters — each thread writes to its own res index, no race condition.
    # Counts how many subsamples hit each boundary for that gene.
    upperBoundaryHits = zeros(Int, totResponses)
    lowerBoundaryHits = zeros(Int, totResponses)

    @info "bEBIC: estimating per-subsample EBIC-optimal lambdas ($totSS subsamples per gene)"
    Threads.@threads for res in ProgressBar(1:totResponses)
        predInds      = responsePredInds[res]
        currPredNum   = length(predInds)
        penaltyFactor = grnData.penaltyMat[res, predInds]
        p             = currPredNum

        bestLambdasSS = zeros(totSS)
        ssSelections  = zeros(totSS, currPredNum)   # binary selection at each subsample's EBIC-optimal lambda
        upperHits     = 0
        lowerHits     = 0

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
                # Log-spaced alternative — swap in if linear grid undersamples the optimal region:
                # lambdaGrid = reverse(exp.(range(log(minLambda), log(maxLambda), length = totLambdas)))
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

            # Boundary detection (only meaningful on a user-supplied grid)
            if minLambda !== nothing && maxLambda !== nothing
                if bestInd == 1
                    upperHits += 1   # optimum at most-regularised end → maxLambda too small
                elseif bestInd == length(lsoln.lambda)
                    lowerHits += 1   # optimum at least-regularised end → minLambda too large
                end
            end
        end

        upperBoundaryHits[res] = upperHits
        lowerBoundaryHits[res] = lowerHits

        # Selection count = number of subsamples in which each edge was selected.
        # Stored as a raw integer count (0–totSS), matching bStARS's stabilityMat
        # convention so that confidence scores are on the same scale across methods.
        selectionCount = vec(sum(ssSelections, dims=1))     # count  (0–totSS)

        networkStability[res, predInds] = selectionCount
        lambdaSSMat[res, :]             = bestLambdasSS
    end

    # Boundary warnings — fired once after the loop with aggregate counts
    if minLambda !== nothing && maxLambda !== nothing
        totUpperFits = sum(upperBoundaryHits)
        totLowerFits = sum(lowerBoundaryHits)
        nUpperGenes  = sum(upperBoundaryHits .> 0)
        nLowerGenes  = sum(lowerBoundaryHits .> 0)
        if totUpperFits > 0
            @warn "bebicEstimateLambdas!: EBIC minimum at upper boundary (maxLambda=$maxLambda) " *
                  "in $totUpperFits gene×subsample fits ($nUpperGenes / $totResponses genes affected). " *
                  "True optimum may lie above maxLambda — consider increasing maxLambda or " *
                  "setting it to nothing for automatic path selection."
        end
        if totLowerFits > 0
            @warn "bebicEstimateLambdas!: EBIC minimum at lower boundary (minLambda=$minLambda) " *
                  "in $totLowerFits gene×subsample fits ($nLowerGenes / $totResponses genes affected). " *
                  "True optimum may lie below minLambda — consider decreasing minLambda or " *
                  "setting it to nothing for automatic path selection."
        end
    end

    # Zero-out Inf entries
    nsVec                = networkStability[:]
    nsVec[isinf.(nsVec)] .= 0.0
    networkStability     = reshape(nsVec, totResponses, totPredictors)

    buildGrn.networkStability = networkStability
    grnData.lambdaSS          = lambdaSSMat

    @info "bEBIC: subsample lambda estimation complete" genesProcessed=totResponses totSubsamples=totSS
end


"""
    bebicFinalFit!(grnData, buildGrn; kwargs...)

Aggregate the per-subsample EBIC-optimal lambdas (stored in `grnData.lambdaSS`)
into one lambda per gene, then perform a final full-data LASSO fit at that lambda.

Must be called after `bebicEstimateLambdas!`.

# Arguments
- `grnData`     : `GrnData` with `lambdaSS` populated by `bebicEstimateLambdas!`
- `buildGrn`    : `BuildGrn` with `networkStability` already filled
- `aggregateFn` : Function applied row-wise to `lambdaSS` to produce one lambda
                  per gene (default `median`). Pass any `AbstractVector → Real`
                  function, e.g. `mean`, `x -> quantile(x, 0.25)`.
- `zScoreLASSO` : Z-score predictors before fitting (default true; should match
                  the value used in `bebicEstimateLambdas!`)
- `outputDir`   : If set, saves `bebicOutMat.jld` (full `GrnData` struct, preserves
                  `lambdaSS`) and `bebic_lambda_summary.tsv` (per-gene median ± std).
                  Default `nothing`.

# Fills
- `buildGrn.betas`      : signed coefficients at the aggregated lambda per gene
- `buildGrn.lambda`     : aggregated lambda per gene (one value per response)
- `grnData.ebicLambdas` : same as `buildGrn.lambda`

# Notes
- Only the aggregation strategy (`aggregateFn`) and the final fit are performed
  here. Re-running with a different `aggregateFn` (e.g. `mean` instead of
  `median`) is cheap — the expensive subsample fits in `bebicEstimateLambdas!`
  do not need to be repeated.
- `gamma` does not appear here because EBIC scoring (which uses `gamma`) is done
  in `bebicEstimateLambdas!`.
"""
function bebicFinalFit!(
    grnData::GrnData,
    buildGrn::BuildGrn;
    aggregateFn::Function             = median,
    zScoreLASSO::Bool                 = true,
    outputDir::Union{String, Nothing} = nothing
)
    totResponses  = size(grnData.responseMat, 1)
    totPredictors = size(grnData.predictorMat, 1)

    if isempty(grnData.lambdaSS)
        error("lambdaSS is empty. Call bebicEstimateLambdas! before bebicFinalFit!.")
    end

    responsePredInds = [findall(x -> x != Inf, grnData.penaltyMat[r, :])
                        for r in 1:totResponses]

    # Aggregate per-subsample lambdas → one per gene
    ebicLambdas = vec(mapslices(aggregateFn, grnData.lambdaSS, dims=2))
    betas_out   = zeros(totResponses, totPredictors)

    Threads.@threads for res in 1:totResponses
        predInds      = responsePredInds[res]
        penaltyFactor = grnData.penaltyMat[res, predInds]
        finalLambda   = ebicLambdas[res]

        dtFull    = fit(ZScoreTransform, grnData.predictorMat[predInds, :], dims=2)
        fullPreds = transpose(StatsBase.transform(dtFull, grnData.predictorMat[predInds, :]))

        if zScoreLASSO
            dtRFull       = fit(ZScoreTransform, grnData.responseMat[res, :], dims=1)
            fullResponses = StatsBase.transform(dtRFull, grnData.responseMat[res, :])
        else
            fullResponses = grnData.responseMat[res, :]
        end

        lsoln = glmnet(fullPreds, vec(fullResponses);
                       penalty_factor = penaltyFactor,
                       lambda         = [finalLambda],
                       alpha          = 1.0)
        betas_out[res, predInds] = vec(lsoln.betas[:, 1])
    end

    buildGrn.lambda     = ebicLambdas
    buildGrn.betas      = betas_out
    grnData.ebicLambdas = ebicLambdas

    if !isnothing(outputDir) && outputDir != ""
        # Save full grnData — preserves lambdaSS (nGenes × totSS per-subsample
        # lambda matrix) which cannot be reconstructed from the TSV summary alone.
        save_object(joinpath(outputDir, "bebicOutMat.jld"), grnData)
        # Write per-gene lambda summary TSV (median + std across subsamples).
        writeBEBICLambdaTable!(grnData, buildGrn; outputDir = outputDir)
    end

    totSS = size(grnData.lambdaSS, 2)
    @info "bEBIC complete" genesProcessed=totResponses lambdaMedian=median(ebicLambdas) lambdaMin=minimum(ebicLambdas) lambdaMax=maximum(ebicLambdas) totSubsamples=totSS
end
