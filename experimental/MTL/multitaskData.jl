"""
multitaskData.jl  [NEW FILE]

Defines the data structures needed for multitask inference.
The key addition over the original codebase is splitting one expression
matrix into per-task views, and storing per-task GRN results.
"""

# ── Task-split expression data ────────────────────────────────────────────────
"""
    MultitaskExpressionData

Holds per-task views of the expression data. Each task corresponds to
a biological context (e.g. cell type) defined by the column labels of
the original expression matrix.

Fields
──────
- `tasks`          : names of each task (e.g. ["Tfh", "Th1", "Treg"])
- `taskData`       : one GeneExpressionData per task            [per-task]
- `taskSimilarity` : tasks × tasks similarity matrix            [NEW]
- `taskLabels`     : maps each sample column → task name        [NEW]
"""
mutable struct MultitaskExpressionData
    tasks::Vector{String}
    taskData::Vector{GeneExpressionData}
    taskSimilarity::Matrix{Float64}
    taskLabels::Vector{String}

    function MultitaskExpressionData()
        return new(
            String[],
            GeneExpressionData[],
            Matrix{Float64}(undef, 0, 0),
            String[]
        )
    end
end


# ── Per-task GRN data ─────────────────────────────────────────────────────────
"""
    MultitaskGrnData

Wraps one GrnData per task plus the fusion penalty parameter.
The per-task GrnData structs are identical to the original codebase —
the only addition is `fusionLambda` which controls how strongly
related tasks are pulled toward each other.

Fields
──────
- `taskGrnData`   : one GrnData per task                        [per-task]
- `fusionLambda`  : λ_f — fusion penalty strength               [NEW]
- `taskGraph`     : tasks × tasks similarity (copied from MT data for convenience)
"""
mutable struct MultitaskGrnData
    taskGrnData::Vector{GrnData}
    fusionLambda::Float64
    taskGraph::Matrix{Float64}

    function MultitaskGrnData()
        return new(
            GrnData[],
            0.1,
            Matrix{Float64}(undef, 0, 0)
        )
    end
end


# ── Per-task BuildGrn ─────────────────────────────────────────────────────────
"""
    MultitaskBuildGrn

Holds one BuildGrn per task plus the consensus network 
(averaged across tasks after fusion).

Fields
──────
- `taskBuildGrn`      : one BuildGrn per task                   [per-task]
- `consensusNetwork`  : gene × TF matrix averaged across tasks  [NEW]
- `tasks`             : task names for indexing                 [NEW]
"""
mutable struct MultitaskBuildGrn
    taskBuildGrn::Vector{BuildGrn}
    consensusNetwork::Matrix{Float64}
    tasks::Vector{String}

    function MultitaskBuildGrn()
        return new(
            BuildGrn[],
            Matrix{Float64}(undef, 0, 0),
            String[]
        )
    end
end
