"""
taskSimilarity.jl  [NEW FILE]

Constructs the task similarity graph used as the fusion penalty structure.
The graph determines which tasks are pulled toward each other by λ_f.

Key design principle: similarity should come from an INDEPENDENT source,
not from the same expression data used to fit the model (circular).
Three options are provided in order of preference.
"""


"""
    buildSimilarityFromMetadata(tasks, metadataFile; connector="_")

Build task similarity from a user-supplied metadata file.
This is the preferred approach — avoids circular use of expression data.

The metadata file should be a TSV with columns: Task, Group
where Group defines which tasks are biologically related.
Tasks in the same group get similarity = 1, across groups = 0.

# Arguments
- `tasks`        : task names in order matching MultitaskExpressionData.tasks
- `metadataFile` : path to TSV metadata file
- `connector`    : separator used in task names (default "_")
"""
function buildSimilarityFromMetadata(tasks::Vector{String}, metadataFile::String)
    df = CSV.read(metadataFile, DataFrame; delim='\t')
    nTasks = length(tasks)
    S = zeros(Float64, nTasks, nTasks)

    taskToGroup = Dict(row.Task => row.Group for row in eachrow(df))

    for i in 1:nTasks
        for j in 1:nTasks
            gi = get(taskToGroup, tasks[i], nothing)
            gj = get(taskToGroup, tasks[j], nothing)
            if gi !== nothing && gj !== nothing && gi == gj
                S[i, j] = 1.0
            end
        end
    end
    # diagonal is always 1
    for i in 1:nTasks
        S[i, i] = 1.0
    end
    return S
end


"""
    buildSimilarityFromExpression(mtData::MultitaskExpressionData; method=:pearson)

Build task similarity from mean expression profiles per task.

⚠️  WARNING: This is circular — you are using the same expression data
to define task similarity AND to fit the model. Only use this as a 
last resort when no independent metadata is available. Consider 
using only a held-out subset of genes (e.g. housekeeping genes) 
for the similarity calculation.

# Arguments
- `mtData`  : MultitaskExpressionData with taskData populated
- `method`  : :pearson (default) or :spearman
"""
function buildSimilarityFromExpression(mtData::MultitaskExpressionData; method::Symbol=:pearson)
    nTasks = length(mtData.tasks)
    # compute mean expression profile per task
    meanProfiles = hcat([vec(mean(td.targGeneMat, dims=2)) for td in mtData.taskData]...)  # genes × tasks
    S = zeros(Float64, nTasks, nTasks)

    for i in 1:nTasks
        for j in 1:nTasks
            if method == :pearson
                S[i, j] = cor(meanProfiles[:, i], meanProfiles[:, j])
            elseif method == :spearman
                ri = tiedrank(meanProfiles[:, i])
                rj = tiedrank(meanProfiles[:, j])
                S[i, j] = cor(ri, rj)
            end
        end
    end

    # clip to [0, 1] — negative correlations treated as no similarity
    S = max.(S, 0.0)
    return S
end


"""
    buildSimilarityFromOntology(tasks, ontologyFile)

Build task similarity from a cell type ontology distance file.
The ontology file should be a TSV with columns: Task1, Task2, Distance
where Distance is a non-negative value (0 = identical, larger = more distant).

Similarity is computed as sim = exp(-distance / scale) where scale
is the median pairwise distance.
"""
function buildSimilarityFromOntology(tasks::Vector{String}, ontologyFile::String)
    df = CSV.read(ontologyFile, DataFrame; delim='\t')
    nTasks = length(tasks)
    distMat = zeros(Float64, nTasks, nTasks)
    taskIdx = Dict(t => i for (i, t) in enumerate(tasks))

    for row in eachrow(df)
        if haskey(taskIdx, row.Task1) && haskey(taskIdx, row.Task2)
            i, j = taskIdx[row.Task1], taskIdx[row.Task2]
            distMat[i, j] = row.Distance
            distMat[j, i] = row.Distance
        end
    end

    # convert distance to similarity via RBF kernel
    offDiag = [distMat[i,j] for i in 1:nTasks for j in 1:nTasks if i != j]
    scale = isempty(offDiag) ? 1.0 : median(offDiag)
    S = exp.(-distMat ./ max(scale, 1e-8))
    return S
end


"""
    normalizeSimilarity!(S::Matrix{Float64})

Row-normalize similarity matrix so rows sum to 1.
This ensures the fusion penalty is scale-invariant across tasks
with different numbers of neighbors.
"""
function normalizeSimilarity!(S::Matrix{Float64})
    for i in 1:size(S, 1)
        rowSum = sum(S[i, :])
        if rowSum > 0
            S[i, :] ./= rowSum
        end
    end
    return S
end
