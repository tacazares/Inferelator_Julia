"""
admm.jl  [UPDATED]

Changes vs original:
─────────────────────────────────────────────────────────────────────
NEW  elasticNetAlpha parameter throughout (default 1.0 = pure LASSO)
NEW  selectLambdas()           — top-level dispatcher for all three strategies
NEW  selectLambdas_fixedRatio  — Option 1: lambda_f = fusionRatio * lambda_s
NEW  selectLambdas_ebic        — Option 2: EBIC grid search
NEW  selectLambdas_bstars2d    — Option 3: 2D bStARS instability surface
"""

@inline function softThreshold(x::Float64, threshold::Float64)
    return sign(x) * max(abs(x) - threshold, 0.0)
end


# ────────────────────────────────────────────────────────────────────────────
# Core ADMM solver  [UPDATED: elasticNetAlpha added]
# ────────────────────────────────────────────────────────────────────────────

"""
    admm_fused_lasso(predictorMats, responseMats, penaltyFactors,
                     taskSimilarity, lambda_s, lambda_f;
                     rho, maxIter, tol, elasticNetAlpha)

Solve the graph-fused multitask LASSO for a single target gene.

elasticNetAlpha: L1/L2 mixing parameter passed to GLMNet alpha argument.
  1.0 (default) = pure LASSO
  0.5           = equal L1+L2 — recommended when demerged TFs are present
  0.0           = pure ridge
"""
function admm_fused_lasso(
        predictorMats::Vector{Matrix{Float64}},
        responseMats::Vector{Vector{Float64}},
        penaltyFactors::Vector{Vector{Float64}},
        taskSimilarity::Matrix{Float64},
        lambda_s::Float64,
        lambda_f::Float64;
        rho::Float64             = 1.0,
        maxIter::Int             = 100,
        tol::Float64             = 1e-4,
        elasticNetAlpha::Float64 = 1.0
    )

    nTasks = length(predictorMats)
    nTFs   = size(predictorMats[1], 1)
    edges  = [(d, dp) for d in 1:nTasks for dp in (d+1):nTasks
              if taskSimilarity[d, dp] > 0]
    nEdges = length(edges)

    W = zeros(Float64, nTFs, nTasks)
    Z = zeros(Float64, nTFs, nEdges)
    U = zeros(Float64, nTFs, nEdges)

    for iter in 1:maxIter
        W_prev = copy(W)

        # W-update: per-task elastic net with ADMM augmentation
        for d in 1:nTasks
            dEdges  = [(e, d1 == d ? 1.0 : -1.0)
                       for (e, (d1, d2)) in enumerate(edges)
                       if d1 == d || d2 == d]
            nSamps  = size(predictorMats[d], 2)
            A       = transpose(predictorMats[d])
            x       = responseMats[d]

            if isempty(dEdges)
                lsoln = glmnet(A, x,
                               penalty_factor = penaltyFactors[d],
                               lambda         = [lambda_s],
                               alpha          = elasticNetAlpha)
                W[:, d] = vec(lsoln.betas)
            else
                augA = copy(A)
                augX = copy(x)
                for (e, sgn) in dEdges
                    target  = sgn > 0 ? Z[:, e] - U[:, e] : -Z[:, e] - U[:, e]
                    sqrtRho = sqrt(rho)
                    augA    = vcat(augA, sqrtRho * I(nTFs))
                    augX    = vcat(augX, sqrtRho * target)
                end
                lsoln = glmnet(augA, augX,
                               penalty_factor = penaltyFactors[d],
                               lambda         = [lambda_s / (nSamps + nTFs * length(dEdges))],
                               alpha          = elasticNetAlpha)
                W[:, d] = vec(lsoln.betas)
            end
        end

        # Z-update: fusion proximal operator
        for (e, (d, dp)) in enumerate(edges)
            diff      = W[:, d] - W[:, dp] + U[:, e]
            threshold = lambda_f * taskSimilarity[d, dp] / rho
            Z[:, e]   = softThreshold.(diff, threshold)
        end

        # U-update: dual variable
        for (e, (d, dp)) in enumerate(edges)
            U[:, e] += W[:, d] - W[:, dp] - Z[:, e]
        end

        norm(W - W_prev) < tol && break
    end

    return W
end


function admmWarmStart(
        predictorMats::Vector{Matrix{Float64}},
        responseMats::Vector{Vector{Float64}},
        penaltyFactors::Vector{Vector{Float64}},
        taskSimilarity::Matrix{Float64},
        lambdaRange_s::Vector{Float64},
        lambda_f::Float64;
        elasticNetAlpha::Float64 = 1.0,
        kwargs...
    )
    nTFs    = size(predictorMats[1], 1)
    nTasks  = length(predictorMats)
    nLambda = length(lambdaRange_s)
    betasByLambda = Array{Float64, 3}(undef, nTFs, nTasks, nLambda)

    for (li, ls) in enumerate(lambdaRange_s)
        betasByLambda[:, :, li] = admm_fused_lasso(
            predictorMats, responseMats, penaltyFactors,
            taskSimilarity, ls, lambda_f;
            elasticNetAlpha = elasticNetAlpha, kwargs...
        )
    end
    return betasByLambda
end


# ────────────────────────────────────────────────────────────────────────────
# Lambda selection dispatcher  [NEW]
# ────────────────────────────────────────────────────────────────────────────

"""
    selectLambdas(mtGrnData, mtData, res; lambdaOpt, ...)

Top-level dispatcher for joint (lambda_s, lambda_f) selection.

lambdaOpt:
  :fixed_ratio  — lambda_f = fusionRatio * lambda_s  (default, zero added cost)
  :ebic         — EBIC grid search over 2D (lambda_s, lambda_f)
  :bstars_2d    — 2D bStARS instability surface + EBIC tiebreaker

Returns (lambda_s, lambda_f) for target gene `res`.
"""
function selectLambdas(
        mtGrnData::MultitaskGrnData,
        mtData::MultitaskExpressionData,
        res::Int;
        lambdaOpt::Symbol          = :fixed_ratio,
        fusionRatio::Float64       = 0.1,
        ebicGamma::Float64         = 1.0,
        gridSize::Int              = 10,
        refinementSize::Int        = 10,
        targetInstability::Float64 = 0.05,
        elasticNetAlpha::Float64   = 1.0,
        totSS::Int                 = 20,
        zTarget::Bool              = false
    )

    if lambdaOpt == :fixed_ratio
        return _selectLambdas_fixedRatio(mtGrnData, res; fusionRatio)
    elseif lambdaOpt == :ebic
        return _selectLambdas_ebic(mtGrnData, mtData, res;
                                    ebicGamma, gridSize, elasticNetAlpha, totSS, zTarget)
    elseif lambdaOpt == :bstars_2d
        return _selectLambdas_bstars2d(mtGrnData, mtData, res;
                                        targetInstability, ebicGamma,
                                        coarseSize = gridSize,
                                        fineSize   = refinementSize,
                                        elasticNetAlpha, totSS, zTarget)
    else
        error("Unknown lambdaOpt: $lambdaOpt. Choose :fixed_ratio, :ebic, or :bstars_2d")
    end
end


# ── Option 1: Fixed ratio ────────────────────────────────────────────────────
function _selectLambdas_fixedRatio(mtGrnData::MultitaskGrnData, res::Int;
                                    fusionRatio::Float64 = 0.1)
    refGrn          = mtGrnData.taskGrnData[1]
    lambdaRangeGene = refGrn.lambdaRangeGene[res]
    currInstabs     = refGrn.geneInstabilities[res, :]
    devs            = abs.(currInstabs .- 0.05)
    globalMin       = minimum(devs)
    minInd          = findall(x -> x == globalMin, devs)[end]
    lambda_s        = lambdaRangeGene[minInd]
    lambda_f        = fusionRatio * lambda_s
    return lambda_s, lambda_f
end


# ── Option 2: EBIC grid search ───────────────────────────────────────────────
function _selectLambdas_ebic(
        mtGrnData::MultitaskGrnData,
        mtData::MultitaskExpressionData,
        res::Int;
        ebicGamma::Float64       = 1.0,
        gridSize::Int            = 10,
        elasticNetAlpha::Float64 = 1.0,
        totSS::Int               = 20,
        zTarget::Bool            = false
    )
    nTasks           = length(mtData.tasks)
    responsePredInds = [findall(x -> x != Inf, mtGrnData.taskGrnData[d].penaltyMat[res, :])
                        for d in 1:nTasks]
    taskPredMats     = [transpose(mtGrnData.taskGrnData[d].predictorMat[responsePredInds[d], :])
                        for d in 1:nTasks]
    taskRespVecs     = [vec(mtGrnData.taskGrnData[d].responseMat[res, :]) for d in 1:nTasks]
    taskPenalties    = [mtGrnData.taskGrnData[d].penaltyMat[res, responsePredInds[d]]
                        for d in 1:nTasks]

    lambda_s_grid = exp.(range(log(0.01), log(10.0), length = gridSize))
    alpha_grid    = range(0.05, 1.0, length = gridSize)

    # Traverse from high to low (lambda_s + lambda_f) for warm starting
    grid_pairs = [(ls, a * ls)
                  for ls in reverse(lambda_s_grid)
                  for a  in reverse(alpha_grid)
                  if 0.5 < ls / (a * ls) < 2.0]

    bestEBIC     = Inf
    best_ls      = grid_pairs[1][1]
    best_lf      = grid_pairs[1][2]

    for (ls, lf) in grid_pairs
        W    = admm_fused_lasso(taskPredMats, taskRespVecs, taskPenalties,
                                 mtGrnData.taskGraph, ls, lf;
                                 elasticNetAlpha = elasticNetAlpha)
        ebic = _ebic(W, taskPredMats, taskRespVecs, nTasks, ebicGamma)
        if ebic < bestEBIC
            bestEBIC = ebic; best_ls = ls; best_lf = lf
        end
    end
    return best_ls, best_lf
end


# ── Option 3: 2D bStARS ─────────────────────────────────────────────────────
function _selectLambdas_bstars2d(
        mtGrnData::MultitaskGrnData,
        mtData::MultitaskExpressionData,
        res::Int;
        targetInstability::Float64 = 0.05,
        ebicGamma::Float64         = 1.0,
        coarseSize::Int            = 5,
        fineSize::Int              = 10,
        elasticNetAlpha::Float64   = 1.0,
        totSS::Int                 = 20,
        zTarget::Bool              = false
    )
    nTasks           = length(mtData.tasks)
    responsePredInds = [findall(x -> x != Inf, mtGrnData.taskGrnData[d].penaltyMat[res, :])
                        for d in 1:nTasks]

    # Stage 1: coarse grid
    ls_c = exp.(range(log(0.01), log(5.0),  length = coarseSize))
    lf_c = exp.(range(log(0.001), log(2.0), length = coarseSize))
    D_c  = _instabGrid(mtGrnData, mtData, res, ls_c, lf_c,
                        responsePredInds, nTasks, totSS, zTarget, elasticNetAlpha)

    # Identify region near contour
    near   = findall(x -> abs(x - targetInstability) < targetInstability * 0.5, D_c)
    if isempty(near)
        _, idx = findmin(abs.(D_c .- targetInstability))
        near   = [idx]
    end
    rs = first.(Tuple.(near)); cs = last.(Tuple.(near))
    ls_range = ls_c[max(1, minimum(rs)-1):min(coarseSize, maximum(rs)+1)]
    lf_range = lf_c[max(1, minimum(cs)-1):min(coarseSize, maximum(cs)+1)]

    # Stage 2: fine grid in region
    ls_f = collect(range(minimum(ls_range), maximum(ls_range), length = fineSize))
    lf_f = collect(range(minimum(lf_range), maximum(lf_range), length = fineSize))
    D_f  = _instabGrid(mtGrnData, mtData, res, ls_f, lf_f,
                        responsePredInds, nTasks, totSS, zTarget, elasticNetAlpha)

    # Points on contour
    onContour = findall(x -> abs(x - targetInstability) < targetInstability * 0.3, D_f)
    if isempty(onContour)
        _, idx    = findmin(abs.(D_f .- targetInstability))
        onContour = [idx]
    end

    # EBIC tiebreaker
    taskPredMats  = [transpose(mtGrnData.taskGrnData[d].predictorMat[responsePredInds[d], :])
                     for d in 1:nTasks]
    taskRespVecs  = [vec(mtGrnData.taskGrnData[d].responseMat[res, :]) for d in 1:nTasks]
    taskPenalties = [mtGrnData.taskGrnData[d].penaltyMat[res, responsePredInds[d]]
                     for d in 1:nTasks]

    bestEBIC = Inf; best_ls = ls_f[1]; best_lf = lf_f[1]
    for idx in onContour
        r, c = Tuple(idx)
        ls   = ls_f[r]; lf = lf_f[c]
        W    = admm_fused_lasso(taskPredMats, taskRespVecs, taskPenalties,
                                 mtGrnData.taskGraph, ls, lf;
                                 elasticNetAlpha = elasticNetAlpha)
        ebic = _ebic(W, taskPredMats, taskRespVecs, nTasks, ebicGamma)
        if ebic < bestEBIC
            bestEBIC = ebic; best_ls = ls; best_lf = lf
        end
    end
    return best_ls, best_lf
end


# ── Shared helpers ───────────────────────────────────────────────────────────
function _ebic(W::Matrix{Float64},
               taskPredMats::Vector{Matrix{Float64}},
               taskRespVecs::Vector{Vector{Float64}},
               nTasks::Int,
               gamma::Float64)
    ebic = 0.0
    for d in 1:nTasks
        n_d  = length(taskRespVecs[d])
        p_d  = size(taskPredMats[d], 2)
        w_d  = W[:, d]
        k_d  = sum(abs.(w_d) .> 1e-10)
        yhat = taskPredMats[d] * w_d
        rss  = max(sum((taskRespVecs[d] .- yhat).^2), 1e-12)
        logC = (k_d > 0 && p_d > k_d) ?
               lgamma(p_d+1) - lgamma(k_d+1) - lgamma(p_d-k_d+1) : 0.0
        ebic += n_d * log(rss / n_d) + k_d * log(n_d) + 2 * gamma * logC
    end
    return ebic / nTasks
end


function _instabGrid(
        mtGrnData, mtData, res,
        ls_grid, lf_grid,
        responsePredInds, nTasks, totSS, zTarget, elasticNetAlpha
    )
    nLS     = length(ls_grid)
    nLF     = length(lf_grid)
    instabs = zeros(Float64, nLS, nLF)
    subsamps = mtGrnData.taskGrnData[1].subsamps

    for (li, ls) in enumerate(ls_grid), (fi, lf) in enumerate(lf_grid)
        ssVals = zeros(Float64, nTasks, length(responsePredInds[1]))
        nss    = min(totSS, size(subsamps, 1))

        for ss in 1:nss
            preds = Matrix{Float64}[]; resps = Vector{Float64}[]; pens = Vector{Float64}[]
            for d in 1:nTasks
                sub      = mtGrnData.taskGrnData[d].subsamps[ss, :]
                pidx     = responsePredInds[d]
                dt       = fit(ZScoreTransform,
                               mtGrnData.taskGrnData[d].predictorMat[pidx, sub], dims=2)
                cp       = transpose(StatsBase.transform(dt,
                               mtGrnData.taskGrnData[d].predictorMat[pidx, sub]))
                cr       = zTarget ?
                    StatsBase.transform(fit(ZScoreTransform,
                        mtGrnData.taskGrnData[d].responseMat[res, sub], dims=1),
                        mtGrnData.taskGrnData[d].responseMat[res, sub]) :
                    mtGrnData.taskGrnData[d].responseMat[res, sub]
                push!(preds, transpose(cp)); push!(resps, vec(cr))
                push!(pens, mtGrnData.taskGrnData[d].penaltyMat[res, pidx])
            end
            W = admm_fused_lasso(preds, resps, pens, mtGrnData.taskGraph, ls, lf;
                                  elasticNetAlpha = elasticNetAlpha)
            for d in 1:nTasks
                ssVals[d, :] += vec(sum(abs.(sign.(W)), dims=2))
            end
        end

        theta2        = ssVals ./ nss
        instabPerEdge = 2 .* theta2 .* (1 .- theta2)
        instabs[li, fi] = mean(instabPerEdge)
    end
    return instabs
end
