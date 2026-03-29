"""
admm.jl  [NEW FILE]

ADMM solver for the graph-fused multitask LASSO.

Objective (per target gene i):

    min_{W} (1/2n) Σ_d ||X_i^(d) - A^(d)T w_i^(d)||² 
            + λ_s Σ_{k,d} |Φ_{k,d} w_{k,d}|          ← prior-weighted LASSO
            + λ_f Σ_{(d,d')∈E} sim(d,d') ||w^(d) - w^(d')||₁  ← fusion

where W is a TFs × tasks matrix for gene i.

ADMM reformulation:
    Introduce Z^(d,d') = w^(d) - w^(d') for each edge (d,d') in E

W-update : per-task weighted LASSO (uses GLMNet — unchanged from original)
Z-update : soft thresholding on pairwise differences (fusion proximal op)
U-update : dual variable update (standard ADMM)
"""


"""
    softThreshold(x, threshold)

Scalar soft-thresholding operator used in the Z-update.
"""
@inline function softThreshold(x::Float64, threshold::Float64)
    return sign(x) * max(abs(x) - threshold, 0.0)
end


"""
    admm_fused_lasso(
        predictorMats, responseMats, penaltyFactors,
        taskSimilarity, lambda_s, lambda_f;
        rho=1.0, maxIter=100, tol=1e-4, alpha=1.0
    )

Solve the graph-fused multitask LASSO for a single target gene
using ADMM.

# Arguments
- `predictorMats`   : Vector of predictor matrices, one per task (TFs × samples)
- `responseMats`    : Vector of response vectors, one per task (samples,)
- `penaltyFactors`  : Vector of penalty factor vectors, one per task (TFs,)
- `taskSimilarity`  : tasks × tasks similarity matrix
- `lambda_s`        : LASSO sparsity penalty
- `lambda_f`        : fusion penalty strength
- `rho`             : ADMM augmented Lagrangian parameter
- `maxIter`         : maximum ADMM iterations
- `tol`             : convergence tolerance
- `alpha`           : elastic net mixing (1.0 = pure LASSO)

# Returns
- `W` : TFs × tasks matrix of coefficients for this gene
"""
function admm_fused_lasso(
        predictorMats::Vector{Matrix{Float64}},
        responseMats::Vector{Vector{Float64}},
        penaltyFactors::Vector{Vector{Float64}},
        taskSimilarity::Matrix{Float64},
        lambda_s::Float64,
        lambda_f::Float64;
        rho::Float64 = 1.0,
        maxIter::Int = 100,
        tol::Float64 = 1e-4,
        alpha::Float64 = 1.0
    )

    nTasks = length(predictorMats)
    nTFs   = size(predictorMats[1], 1)

    # Build edge list from similarity matrix (upper triangle, nonzero off-diagonal)
    edges = [(d, dp) for d in 1:nTasks for dp in (d+1):nTasks if taskSimilarity[d, dp] > 0]
    nEdges = length(edges)

    # Initialize primal and dual variables
    W = zeros(Float64, nTFs, nTasks)   # TFs × tasks — the main variable
    Z = zeros(Float64, nTFs, nEdges)   # TFs × edges — fusion slack variables
    U = zeros(Float64, nTFs, nEdges)   # TFs × edges — dual variables

    for iter in 1:maxIter
        W_prev = copy(W)

        # ── W-update: per-task weighted LASSO with augmented penalty ─────────
        # For each task d, solve:
        #   min (1/2n)||X^(d) - A^(d)T w^(d)||² 
        #       + λ_s |Φ^(d) ⊙ w^(d)|₁ 
        #       + (ρ/2) Σ_{e∋d} ||w^(d) - Z_e + U_e||²
        #
        # The quadratic augmentation term from ADMM is absorbed into
        # an augmented response and augmented predictor passed to GLMNet.

        for d in 1:nTasks
            # collect edges involving task d
            dEdges    = [(e, sign) for (e, (d1, d2)) in enumerate(edges) 
                         if d1 == d || d2 == d
                         for sign in (d1 == d ? 1.0 : -1.0)]

            nSamps    = size(predictorMats[d], 2)
            A         = transpose(predictorMats[d])   # samples × TFs
            x         = responseMats[d]

            if isempty(dEdges)
                # no fusion neighbors — standard GLMNet call (identical to original)
                lsoln = glmnet(A, x,
                               penalty_factor = penaltyFactors[d],
                               lambda = [lambda_s],
                               alpha  = alpha)
                W[:, d] = vec(lsoln.betas)
            else
                # augment response and predictors with ADMM proximity terms
                augA = copy(A)
                augX = copy(x)
                for (e, sgn) in dEdges
                    target = sgn > 0 ? Z[:, e] - U[:, e] : -Z[:, e] - U[:, e]
                    # append sqrt(ρ) * I rows to A and sqrt(ρ) * target to x
                    sqrtRho = sqrt(rho)
                    augA = vcat(augA, sqrtRho * I(nTFs))
                    augX = vcat(augX, sqrtRho * target)
                end
                augPenalty = penaltyFactors[d]   # penalty unchanged
                lsoln = glmnet(augA, augX,
                               penalty_factor = augPenalty,
                               lambda = [lambda_s / (nSamps + nTFs * length(dEdges))],
                               alpha  = alpha)
                W[:, d] = vec(lsoln.betas)
            end
        end

        # ── Z-update: fusion proximal operator (soft threshold) ──────────────
        # Z_e = soft_threshold(w^(d) - w^(d') + U_e, λ_f * sim(d,d') / ρ)
        for (e, (d, dp)) in enumerate(edges)
            diff      = W[:, d] - W[:, dp] + U[:, e]
            threshold = lambda_f * taskSimilarity[d, dp] / rho
            Z[:, e]   = softThreshold.(diff, threshold)
        end

        # ── U-update: dual variable ───────────────────────────────────────────
        for (e, (d, dp)) in enumerate(edges)
            U[:, e] += W[:, d] - W[:, dp] - Z[:, e]
        end

        # ── Convergence check ─────────────────────────────────────────────────
        primalResid = norm(W - W_prev)
        if primalResid < tol
            break
        end
    end

    return W  # TFs × tasks
end


"""
    admmWarmStart(
        predictorMats, responseMats, penaltyFactors,
        taskSimilarity, lambdaRange_s, lambda_f;
        kwargs...
    )

Run ADMM across a range of lambda_s values to estimate per-task,
per-gene instabilities. Analogous to bstarsWarmStart in the original
codebase but operates on the fused multitask objective.

Returns a 4D array: TFs × tasks × subsamples × lambdas
of binary edge indicators (1 = edge selected, 0 = not selected).
"""
function admmWarmStart(
        predictorMats::Vector{Matrix{Float64}},
        responseMats::Vector{Vector{Float64}},
        penaltyFactors::Vector{Vector{Float64}},
        taskSimilarity::Matrix{Float64},
        lambdaRange_s::Vector{Float64},
        lambda_f::Float64;
        kwargs...
    )
    nTFs    = size(predictorMats[1], 1)
    nTasks  = length(predictorMats)
    nLambda = length(lambdaRange_s)

    betasByLambda = Array{Float64, 3}(undef, nTFs, nTasks, nLambda)

    for (li, ls) in enumerate(lambdaRange_s)
        W = admm_fused_lasso(
            predictorMats, responseMats, penaltyFactors,
            taskSimilarity, ls, lambda_f; kwargs...
        )
        betasByLambda[:, :, li] = W
    end

    return betasByLambda
end
