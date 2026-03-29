"""
prepareMultitask.jl  [NEW FILE]

Multitask versions of the preparation functions from prepareGRN.jl.

What stays the same vs original:
─────────────────────────────────────────────────────────────────────
UNCHANGED  preparePredictorMat!    — called once per task
UNCHANGED  preparePenaltyMatrix!   — called once per task
UNCHANGED  constructSubsamples     — called once per task
UNCHANGED  bstarsWarmStart         — called once per task (coarse pass)
NEW        splitByTask!            — splits expression matrix by task
NEW        preparePredictorMatMT!  — loops preparePredictorMat! over tasks
NEW        preparePenaltyMatrixMT! — loops preparePenaltyMatrix! over tasks
NEW        constructSubsamplesMT!  — loops constructSubsamples over tasks
NEW        bstartsEstimateInstabilityMT — replaces bstartsEstimateInstability,
                                          uses ADMM instead of GLMNet per gene
"""


# ── NEW: Split expression matrix by task ─────────────────────────────────────
"""
    splitByTask!(mtData::MultitaskExpressionData, data::GeneExpressionData, taskLabels::Vector{String})

Split one GeneExpressionData into per-task views based on column labels.

`taskLabels` must have the same length as `data.cellLabels` and maps
each sample column to a task name (e.g. "Tfh", "Th1", "Treg").

This is the only entry point that creates per-task data — everything
downstream just iterates over `mtData.taskData`.
"""
function splitByTask!(mtData::MultitaskExpressionData,
                      data::GeneExpressionData,
                      taskLabels::Vector{String})

    if length(taskLabels) != length(data.cellLabels)
        error("taskLabels length ($(length(taskLabels))) must match number of samples ($(length(data.cellLabels)))")
    end

    uniqueTasks = unique(taskLabels)
    mtData.taskLabels = taskLabels
    mtData.tasks = uniqueTasks

    for task in uniqueTasks
        taskInds = findall(x -> x == task, taskLabels)

        td = GeneExpressionData()
        # shared metadata — same across all tasks
        td.geneNames        = data.geneNames
        td.targGenes        = data.targGenes
        td.potRegs          = data.potRegs
        td.potRegsmRNA      = data.potRegsmRNA
        td.tfaGenes         = data.tfaGenes
        # per-task column subsets
        td.cellLabels       = data.cellLabels[taskInds]
        td.geneExpressionMat = data.geneExpressionMat[:, taskInds]
        td.targGeneMat      = data.targGeneMat[:, taskInds]
        td.potRegMatmRNA    = data.potRegMatmRNA[:, taskInds]
        td.tfaGeneMat       = data.tfaGeneMat[:, taskInds]

        push!(mtData.taskData, td)
    end

    println("Split into $(length(uniqueTasks)) tasks: ", join(uniqueTasks, ", "))
end


# ── NEW: Per-task TFA estimation ──────────────────────────────────────────────
"""
    calculateTFAPerTask!(tfaDataVec, mtData; edgeSS=0, zTarget=false, outputDir=nothing)

Run calculateTFA! independently for each task.
Returns a Vector{PriorTFAData}, one per task.

Each task gets its own TFA estimate because TF activity can differ
substantially across cell types — pooling them would obscure this.
"""
function calculateTFAPerTask!(tfaDataVec::Vector{PriorTFAData},
                               mtData::MultitaskExpressionData;
                               edgeSS::Int = 0,
                               zTarget::Bool = false,
                               outputDir::Union{String, Nothing} = nothing)

    for (d, taskName) in enumerate(mtData.tasks)
        taskDir = outputDir !== nothing ? joinpath(outputDir, taskName) : nothing
        if taskDir !== nothing
            mkpath(taskDir)
        end
        calculateTFA!(tfaDataVec[d], mtData.taskData[d];
                      edgeSS = edgeSS, zscoreTargExp = zTarget, outputDir = taskDir)
        println("TFA estimated for task: ", taskName)
    end
end


# ── NEW: Loop prepare functions over tasks ────────────────────────────────────
"""
    preparePredictorMatMT!(mtGrnData, mtData, tfaDataVec, tfaOpt)

Call preparePredictorMat! for each task. [UNCHANGED per-task logic]
"""
function preparePredictorMatMT!(mtGrnData::MultitaskGrnData,
                                 mtData::MultitaskExpressionData,
                                 tfaDataVec::Vector{PriorTFAData},
                                 tfaOpt::String)
    for d in 1:length(mtData.tasks)
        preparePredictorMat!(mtGrnData.taskGrnData[d], mtData.taskData[d], tfaDataVec[d], tfaOpt)
    end
end


"""
    preparePenaltyMatrixMT!(mtData, mtGrnData, priorFilePenalties, lambdaBias, tfaOpt)

Call preparePenaltyMatrix! for each task. [UNCHANGED per-task logic]
"""
function preparePenaltyMatrixMT!(mtData::MultitaskExpressionData,
                                  mtGrnData::MultitaskGrnData,
                                  priorFilePenalties::Vector{String},
                                  lambdaBias::Vector{Float64},
                                  tfaOpt::String)
    for d in 1:length(mtData.tasks)
        preparePenaltyMatrix!(mtData.taskData[d], mtGrnData.taskGrnData[d],
                               priorFilePenalties, lambdaBias, tfaOpt)
    end
end


"""
    constructSubsamplesMT!(mtData, mtGrnData; totSS, subsampleFrac)

Call constructSubsamples for each task. [UNCHANGED per-task logic]
"""
function constructSubsamplesMT!(mtData::MultitaskExpressionData,
                                  mtGrnData::MultitaskGrnData;
                                  totSS::Int = 50,
                                  subsampleFrac::Float64 = 0.68)
    for d in 1:length(mtData.tasks)
        constructSubsamples(mtData.taskData[d], mtGrnData.taskGrnData[d];
                             totSS = totSS, subsampleFrac = subsampleFrac)
    end
end


# ── NEW: Multitask instability estimation (replaces bstartsEstimateInstability) ──
"""
    bstartsEstimateInstabilityMT!(mtGrnData, mtData; 
                                   totLambdas, instabilityLevel, zTarget, outputDir)

Replaces bstartsEstimateInstability from the original codebase.

Key difference: instead of fitting GLMNet independently per task,
this calls admm_fused_lasso which couples related tasks via λ_f.

The per-gene parallelism (Threads.@threads) is preserved — tasks are
coupled WITHIN each gene's ADMM problem, not across genes.

Instability is computed the same way as the original:
    θ = (1/totSS) * edge_selection_count
    instab = 2 * θ * (1 - θ)
But now per task, giving a tasks × lambdas instability surface per gene.
"""
function bstartsEstimateInstabilityMT!(mtGrnData::MultitaskGrnData,
                                        mtData::MultitaskExpressionData;
                                        totLambdas::Int = 10,
                                        instabilityLevel::String = "Gene",
                                        zTarget::Bool = false,
                                        outputDir::Union{String, Nothing} = nothing)

    nTasks = length(mtData.tasks)
    # use first task to get dimensions — all tasks share same gene/TF sets
    refGrn = mtGrnData.taskGrnData[1]
    totResponses, totSamps = size(refGrn.responseMat)
    totPreds = size(refGrn.predictorMat, 1)
    totSS = size(refGrn.subsamps, 1)

    # Lambda range: use network-level bounds from first task as reference
    # TODO: consider per-task lambda ranges for heterogeneous tasks
    minLambda = minimum([g.minLambdaNet for g in mtGrnData.taskGrnData])
    maxLambda = maximum([g.maxLambdaNet for g in mtGrnData.taskGrnData])
    lambdaRange = reverse(collect(range(minLambda, stop=maxLambda, length=totLambdas)))

    # Store edge selection counts: lambdas × genes × TFs × tasks
    ssMatrix = Inf * ones(totLambdas, totResponses, totPreds, nTasks)
    betas    = Array{Float64, 4}(undef, totResponses, totPreds, totLambdas, nTasks)

    # get finite predictor indices per task per response
    responsePredInds = [[findall(x -> x != Inf, mtGrnData.taskGrnData[d].penaltyMat[res, :])
                          for res in 1:totResponses]
                         for d in 1:nTasks]

    Threads.@threads for res in ProgressBar(1:totResponses)
        for ss in 1:totSS
            # collect per-task predictors and responses for this subsample
            taskPredMats  = Matrix{Float64}[]
            taskRespVecs  = Vector{Float64}[]
            taskPenalties = Vector{Float64}[]

            for d in 1:nTasks
                subsamp  = mtGrnData.taskGrnData[d].subsamps[ss, :]
                predInds = responsePredInds[d][res]

                dt = fit(ZScoreTransform,
                         mtGrnData.taskGrnData[d].predictorMat[predInds, subsamp], dims=2)
                currPreds = transpose(StatsBase.transform(dt,
                             mtGrnData.taskGrnData[d].predictorMat[predInds, subsamp]))

                if zTarget
                    dt2 = fit(ZScoreTransform,
                              mtGrnData.taskGrnData[d].responseMat[res, subsamp], dims=1)
                    currResp = StatsBase.transform(dt2,
                                mtGrnData.taskGrnData[d].responseMat[res, subsamp])
                else
                    currResp = mtGrnData.taskGrnData[d].responseMat[res, subsamp]
                end

                push!(taskPredMats,  transpose(currPreds))   # TFs × subsampled_cells
                push!(taskRespVecs,  vec(currResp))
                push!(taskPenalties, mtGrnData.taskGrnData[d].penaltyMat[res, predInds])
            end

            # ── ADMM: couples tasks via fusion penalty ─────────────────────────
            # This is the key difference from the original codebase.
            # Original: glmnet(currPreds, currResponses, ...) per task independently
            # New:      admm_fused_lasso(...) jointly over all tasks
            betasByLambda = admmWarmStart(
                taskPredMats, taskRespVecs, taskPenalties,
                mtGrnData.taskGraph, lambdaRange, mtGrnData.fusionLambda
            )  # TFs × tasks × lambdas

            for d in 1:nTasks
                predInds = responsePredInds[d][res]
                for (li, _) in enumerate(lambdaRange)
                    ssMatrix[li, res, predInds, d] .+= abs.(sign.(betasByLambda[:, d, li]))
                    betas[res, predInds, li, d] = betasByLambda[:, d, li]
                end
            end
        end
    end

    # compute instabilities per task (same formula as original)
    for d in 1:nTasks
        grnD = mtGrnData.taskGrnData[d]
        geneInstabilities = zeros(totResponses, totLambdas)

        for res in 1:totResponses
            predInds    = responsePredInds[d][res]
            ssVals      = ssMatrix[:, res, predInds, d]
            theta2      = (1 / totSS) * ssVals
            instabPerEdge = 2 * (theta2 .* (1 .- theta2))
            aveInstab   = vec(mean(instabPerEdge, dims=2))
            maxUb       = maximum(aveInstab)
            maxUbInd    = findlast(x -> x == maxUb, aveInstab)
            aveInstab[maxUbInd:end] .= maxUb
            geneInstabilities[res, :] = aveInstab
        end

        grnD.geneInstabilities = geneInstabilities
        grnD.lambdaRange       = lambdaRange
        grnD.stabilityMat      = ssMatrix[:, :, :, d]
        grnD.betas             = betas[:, :, :, d]

        if outputDir !== nothing
            taskDir = joinpath(outputDir, mtData.tasks[d])
            mkpath(taskDir)
            save_object(joinpath(taskDir, "instabOutMat.jld"), grnD)
        end
    end

    println("Multitask instability estimation complete.")
end
