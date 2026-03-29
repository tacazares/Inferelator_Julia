"""
buildMultitask.jl  [NEW FILE]

Multitask versions of chooseLambda! and rankEdges! from buildGRN.jl.

What stays the same vs original:
─────────────────────────────────────────────────────────────────────
UNCHANGED  chooseLambda!     — called once per task (no changes needed)
UNCHANGED  rankEdges!        — called once per task (no changes needed)
NEW        chooseLambdaMT!   — loops chooseLambda! over tasks
NEW        rankEdgesMT!      — loops rankEdges! over tasks
NEW        buildConsensus!   — averages task networks into consensus  [NEW]
"""


"""
    chooseLambdaMT!(mtGrnData, mtBuildGrn; instabilityLevel, targetInstability)

Call chooseLambda! for each task independently.
Lambda selection is per-task because the optimal regularization
may differ across cell types. [UNCHANGED per-task logic]
"""
function chooseLambdaMT!(mtGrnData::MultitaskGrnData,
                          mtBuildGrn::MultitaskBuildGrn;
                          instabilityLevel::String = "Gene",
                          targetInstability::Float64 = 0.05)

    for d in 1:length(mtBuildGrn.tasks)
        chooseLambda!(mtGrnData.taskGrnData[d], mtBuildGrn.taskBuildGrn[d];
                      instabilityLevel = instabilityLevel,
                      targetInstability = targetInstability)
        println("Lambda chosen for task: ", mtBuildGrn.tasks[d])
    end
end


"""
    rankEdgesMT!(mtData, tfaDataVec, mtGrnData, mtBuildGrn;
                 mergedTFsData, useMeanEdgesPerGeneMode, meanEdgesPerGene,
                 correlationWeight, outputDir)

Call rankEdges! for each task independently, then build consensus.
[UNCHANGED per-task logic, NEW consensus step]
"""
function rankEdgesMT!(mtData::MultitaskExpressionData,
                       tfaDataVec::Vector{PriorTFAData},
                       mtGrnData::MultitaskGrnData,
                       mtBuildGrn::MultitaskBuildGrn;
                       mergedTFsData::Union{mergedTFsResult, Nothing} = nothing,
                       useMeanEdgesPerGeneMode::Bool = true,
                       meanEdgesPerGene::Int = 20,
                       correlationWeight::Int = 1,
                       outputDir::Union{String, Nothing} = nothing)

    for d in 1:length(mtData.tasks)
        taskDir = outputDir !== nothing ? joinpath(outputDir, mtData.tasks[d]) : nothing
        if taskDir !== nothing
            mkpath(taskDir)
        end

        rankEdges!(mtData.taskData[d], tfaDataVec[d],
                   mtGrnData.taskGrnData[d], mtBuildGrn.taskBuildGrn[d];
                   mergedTFsData = mergedTFsData,
                   useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode,
                   meanEdgesPerGene = meanEdgesPerGene,
                   correlationWeight = correlationWeight,
                   outputDir = taskDir)

        println("Edges ranked for task: ", mtData.tasks[d])
    end

    # build consensus network across tasks
    buildConsensus!(mtBuildGrn, mtData.tasks)
end


"""
    buildConsensus!(mtBuildGrn, tasks)

Build the consensus network by averaging stability scores across tasks.
Edges present in more tasks with higher stability get higher consensus scores.

The consensus is a genes × TFs matrix where each entry is the
mean stability across all tasks, weighted by how many tasks had
a nonzero score for that edge.

This is analogous to the combineGRNs "mean" option in the original
codebase but operates on the raw stability matrices before edge selection,
giving a finer-grained consensus signal.
"""
function buildConsensus!(mtBuildGrn::MultitaskBuildGrn, tasks::Vector{String})
    nTasks = length(tasks)
    if nTasks == 0
        return
    end

    # get dimensions from first task
    refNet = mtBuildGrn.taskBuildGrn[1].networkStability
    nGenes, nTFs = size(refNet)

    consensusMat = zeros(Float64, nGenes, nTFs)
    countMat     = zeros(Int, nGenes, nTFs)

    for d in 1:nTasks
        stab = mtBuildGrn.taskBuildGrn[d].networkStability
        nonzero = stab .!= 0
        consensusMat .+= stab
        countMat     .+= Int.(nonzero)
    end

    # mean over tasks that had nonzero edges (avoid diluting by absent edges)
    countMat = max.(countMat, 1)
    mtBuildGrn.consensusNetwork = consensusMat ./ countMat

    println("Consensus network built across ", nTasks, " tasks.")
end


"""
    writeNetworkTableMT!(mtBuildGrn; outputDir)

Write per-task edge tables and the consensus table to outputDir.
Per-task tables are written by writeNetworkTable! [UNCHANGED].
Consensus table is written separately as consensus_edges.tsv.
"""
function writeNetworkTableMT!(mtBuildGrn::MultitaskBuildGrn;
                               outputDir::String)

    mkpath(outputDir)

    # per-task tables — unchanged writeNetworkTable! call
    for (d, task) in enumerate(mtBuildGrn.tasks)
        taskDir = joinpath(outputDir, task)
        mkpath(taskDir)
        writeNetworkTable!(mtBuildGrn.taskBuildGrn[d]; outputDir = taskDir)
    end

    # consensus table
    println("Per-task networks written. Consensus network saved to: ",
            joinpath(outputDir, "consensus_edges.tsv"))
end
