module MultitaskInferelator

    # ── Unchanged modules from original Inferelator ──────────────────────────
    include("core/geneExpression.jl")       # Data module         [UNCHANGED]
    include("core/mergeDegenerateTFs.jl")   # MergeDegenerate     [UNCHANGED]
    include("core/networkIO.jl")            # NetworkIO           [UNCHANGED]
    include("utils/dataUtils.jl")           # DataUtils           [UNCHANGED]
    include("utils/utilsGRN.jl")            # firstNByGroup       [UNCHANGED]

    # ── Modified modules ──────────────────────────────────────────────────────
    include("core/priorTFA.jl")             # PriorTFA            [MODIFIED - per-task TFA]
    include("core/GRN.jl")                  # GRN structs         [MODIFIED - new structs]

    # ── New multitask modules ─────────────────────────────────────────────────
    include("multitask/multitaskData.jl")   # MT data structs     [NEW]
    include("multitask/taskSimilarity.jl")  # Study graph         [NEW]
    include("multitask/admm.jl")            # ADMM solver         [NEW]
    include("multitask/prepareMultitask.jl")# MT prepare fns      [NEW]
    include("multitask/buildMultitask.jl")  # MT build/rank fns   [NEW]
    include("multitask/combineMultitask.jl")# MT combine fns      [NEW]

    # ── Evaluation (unchanged) ────────────────────────────────────────────────
    include("evaluation/Metric.jl")
    include("evaluation/calcPRinfTRNs.jl")
    include("evaluation/plotSingleUtils.jl")
    include("evaluation/plotBatchMetrics.jl")

end
