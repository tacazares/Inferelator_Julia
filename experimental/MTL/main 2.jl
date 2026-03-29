cd("/data/miraldiNB/Michael/Scripts/GRN/MultitaskInferelator")
using Pkg
Pkg.activate(".")
using Revise
include("src/MultitaskInferelator.jl")
using .MultitaskInferelator

"""
    runMultitaskInferelator(; kwargs...)

Multitask extension of the original Inferelator pipeline.

Changes vs original runInferelator:
────────────────────────────────────────────────────────────────────
UNCHANGED  loadExpressionData!          load expression matrix
UNCHANGED  loadAndFilterTargetGenes!    filter target genes
UNCHANGED  loadPotentialRegulators!     load TF list
UNCHANGED  processTFAGenes!             set TFA gene set
UNCHANGED  mergeDegenerateTFs           merge degenerate TFs
UNCHANGED  processPriorFile!            process prior file
UNCHANGED  preparePenaltyMatrix!        build penalty matrix (per task)
UNCHANGED  constructSubsamples          subsample indices (per task)
UNCHANGED  bstarsWarmStart              coarse lambda bounds (per task)
UNCHANGED  chooseLambda!                lambda selection (per task)
UNCHANGED  rankEdges!                   edge ranking (per task)
UNCHANGED  writeNetworkTable!           write edge tables (per task)
UNCHANGED  combineGRNs                  combine TFA/TFmRNA networks
UNCHANGED  combineGRNS2                 re-estimate TFA on combined

NEW        splitByTask!                 split expression by cell type
NEW        buildSimilarityFrom*         construct task graph
NEW        calculateTFAPerTask!         per-task TFA estimation
NEW        bstartsEstimateInstabilityMT ADMM-fused instability estimation
NEW        buildConsensus!              consensus network across tasks

New parameters vs original:
────────────────────────────────────────────────────────────────────
taskLabelFile   path to TSV mapping samples → task names
fusionLambda    λ_f fusion penalty strength (default 0.1)
similaritySource  :metadata, :expression, or :ontology
similarityFile  path to metadata/ontology file (if needed)
"""
function runMultitaskInferelator(;
    geneExprFile::String,
    targFile::String,
    regFile::String,
    priorFile::String,
    priorFilePenalties::Vector{String},
    taskLabelFile::String,                         # NEW — maps samples → tasks
    tfaGeneFile::String = "",
    outputDir::String,
    tfaOptions::Vector{String} = ["", "TFmRNA"],
    totSS::Int = 80,
    bstarsTotSS::Int = 5,
    subsampleFrac::Float64 = 0.68,
    minLambda::Float64 = 0.01,
    maxLambda::Float64 = 0.5,
    totLambdasBstars::Int = 20,
    totLambdas::Int = 40,
    targetInstability::Float64 = 0.05,
    meanEdgesPerGene::Int = 20,
    correlationWeight::Int = 1,
    minTargets::Int = 3,
    edgeSS::Int = 0,
    lambdaBias::Vector{Float64} = [0.5],
    instabilityLevel::String = "Gene",
    useMeanEdgesPerGeneMode::Bool = true,
    combineOpt::String = "max",
    zTarget::Bool = true,
    fusionLambda::Float64 = 0.1,                  # NEW — fusion penalty (used by :fixed_ratio)
    similaritySource::Symbol = :expression,        # NEW — how to build task graph
    similarityFile::Union{String, Nothing} = nothing, # NEW — metadata/ontology file
    lambdaOpt::Symbol = :fixed_ratio,              # NEW — :fixed_ratio | :ebic | :bstars_2d
    fusionRatio::Float64 = 0.1,                    # NEW — for :fixed_ratio option
    ebicGamma::Float64 = 1.0,                      # NEW — EBIC gamma (use 1.0 for p>>n)
    gridSize::Int = 10,                            # NEW — lambda grid size for :ebic/:bstars_2d
    refinementSize::Int = 10,                      # NEW — fine grid size for :bstars_2d
    elasticNetAlpha::Float64 = 1.0                 # NEW — L1/L2 mix (0.5 recommended with demerged TFs)
)

    # build output directory
    subsamplePct = subsampleFrac * 100
    subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
    lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")
    networkBaseName = "MT_" * lowercase(instabilityLevel) * "Lambda" * lambdaStr * "_" *
                      string(totSS) * "totSS_" * string(meanEdgesPerGene) * "tfsPerGene_subsamplePCT" * subsampleStr
    dirOut = joinpath(outputDir, networkBaseName)
    mkpath(dirOut)

    println("=== Multitask Inferelator Configuration ===")
    println("Output Directory:      ", dirOut)
    println("Expression File:       ", geneExprFile)
    println("Task Label File:       ", taskLabelFile)
    println("Prior File:            ", priorFile)
    println("Fusion Lambda (λ_f):   ", fusionLambda)
    println("Similarity Source:     ", similaritySource)
    println("===========================================")

    # ── STEP 1: Load expression data [UNCHANGED] ──────────────────────────────
    data = GeneExpressionData()
    loadExpressionData!(data, geneExprFile)
    loadAndFilterTargetGenes!(data, targFile; epsilon=0.01)
    loadPotentialRegulators!(data, regFile)
    processTFAGenes!(data, tfaGeneFile; outputDir=dirOut)

    # ── STEP 2: Split by task [NEW] ───────────────────────────────────────────
    taskLabels = vec(readdlm(taskLabelFile, String))  # one label per sample column
    mtData = MultitaskExpressionData()
    splitByTask!(mtData, data, taskLabels)

    # ── STEP 3: Build task similarity graph [NEW] ─────────────────────────────
    if similaritySource == :metadata && similarityFile !== nothing
        S = buildSimilarityFromMetadata(mtData.tasks, similarityFile)
    elseif similaritySource == :ontology && similarityFile !== nothing
        S = buildSimilarityFromOntology(mtData.tasks, similarityFile)
    else
        println("⚠ Using expression-based similarity — see taskSimilarity.jl for caveats.")
        S = buildSimilarityFromExpression(mtData)
    end
    normalizeSimilarity!(S)
    mtData.taskSimilarity = S

    # ── STEP 4: Merge degenerate TFs [UNCHANGED] ──────────────────────────────
    mergedTFsData = mergedTFsResult()
    mergeDegenerateTFs(mergedTFsData, priorFile; fileFormat=2)

    # ── STEP 5: Process prior file [UNCHANGED] ────────────────────────────────
    # Prior is shared across tasks — one prior file, processed once
    tfaDataTemplate = PriorTFAData()
    processPriorFile!(tfaDataTemplate, data, priorFile; mergedTFsData, minTargets=minTargets)

    # ── STEP 6: Per-task TFA estimation [NEW] ────────────────────────────────
    # Each task gets its own TFA estimate
    tfaDataVec = [deepcopy(tfaDataTemplate) for _ in mtData.tasks]
    calculateTFAPerTask!(tfaDataVec, mtData;
                          edgeSS=edgeSS, zTarget=zTarget, outputDir=dirOut)

    # ── STEP 7: Build GRN for each TFA option ────────────────────────────────
    for tfaOpt in tfaOptions
        optName = tfaOpt == "" ? "TFA" : "TFmRNA"
        instabilitiesDir = joinpath(dirOut, optName)
        mkpath(instabilitiesDir)

        # initialize per-task GrnData
        mtGrnData = MultitaskGrnData()
        mtGrnData.fusionLambda = fusionLambda
        mtGrnData.taskGraph = S
        for _ in mtData.tasks
            push!(mtGrnData.taskGrnData, GrnData())
        end

        # prepare matrices — UNCHANGED per-task logic
        preparePredictorMatMT!(mtGrnData, mtData, tfaDataVec, tfaOpt)
        preparePenaltyMatrixMT!(mtData, mtGrnData, priorFilePenalties, lambdaBias, tfaOpt)
        constructSubsamplesMT!(mtData, mtGrnData; totSS=bstarsTotSS, subsampleFrac=subsampleFrac)

        # coarse warm start — run per task independently [UNCHANGED]
        for d in 1:length(mtData.tasks)
            bstarsWarmStart(mtData.taskData[d], tfaDataVec[d], mtGrnData.taskGrnData[d];
                             minLambda=minLambda, maxLambda=maxLambda,
                             totLambdasBstars=totLambdasBstars,
                             targetInstability=targetInstability, zTarget=zTarget)
        end

        constructSubsamplesMT!(mtData, mtGrnData; totSS=totSS, subsampleFrac=subsampleFrac)

        # instability estimation — ADMM-fused [NEW CORE STEP]
        bstartsEstimateInstabilityMT!(mtGrnData, mtData;
                                       totLambdas       = totLambdas,
                                       instabilityLevel = instabilityLevel,
                                       zTarget          = zTarget,
                                       outputDir        = instabilitiesDir,
                                       lambdaOpt        = lambdaOpt,        # NEW
                                       fusionRatio      = fusionRatio,      # NEW
                                       ebicGamma        = ebicGamma,        # NEW
                                       gridSize         = gridSize,         # NEW
                                       refinementSize   = refinementSize,   # NEW
                                       elasticNetAlpha  = elasticNetAlpha)  # NEW

        # lambda selection and edge ranking — UNCHANGED per-task logic
        mtBuildGrn = MultitaskBuildGrn()
        mtBuildGrn.tasks = mtData.tasks
        for _ in mtData.tasks
            push!(mtBuildGrn.taskBuildGrn, BuildGrn())
        end

        chooseLambdaMT!(mtGrnData, mtBuildGrn;
                         instabilityLevel=instabilityLevel,
                         targetInstability=targetInstability)

        rankEdgesMT!(mtData, tfaDataVec, mtGrnData, mtBuildGrn;
                      mergedTFsData=mergedTFsData,
                      useMeanEdgesPerGeneMode=useMeanEdgesPerGeneMode,
                      meanEdgesPerGene=meanEdgesPerGene,
                      correlationWeight=correlationWeight,
                      outputDir=instabilitiesDir)

        writeNetworkTableMT!(mtBuildGrn; outputDir=instabilitiesDir)
    end

    # ── STEP 8: Combine TFA and TFmRNA networks [UNCHANGED] ──────────────────
    # Combine per-task within each TFA option, then across options
    for task in mtData.tasks
        combinedNetDir = joinpath(dirOut, "Combined", task)
        mkpath(combinedNetDir)
        nets2combine = [
            joinpath(dirOut, "TFA",   task, "edges.tsv"),
            joinpath(dirOut, "TFmRNA", task, "edges.tsv")
        ]
        combineGRNs(nets2combine;
                    combineOpt=combineOpt,
                    meanEdgesPerGene=meanEdgesPerGene,
                    useMeanEdgesPerGeneMode=useMeanEdgesPerGeneMode,
                    saveDir=combinedNetDir,
                    saveName=task)

        # re-estimate TFA on combined network [UNCHANGED]
        netsCombinedSparse = joinpath(combinedNetDir, "combined_$(task)_$(combineOpt)_sp.tsv")
        taskIdx = findfirst(x -> x == task, mtData.tasks)
        combineGRNS2(mtData.taskData[taskIdx], mergedTFsData, tfaGeneFile,
                     netsCombinedSparse, edgeSS, minTargets,
                     geneExprFile, targFile, regFile; outputDir=combinedNetDir)
    end

    println("=== Multitask Inferelator Complete ===")
end


# ── Run ────────────────────────────────────────────────────────────────────────
runMultitaskInferelator(
    geneExprFile       = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/counts_Tfh10_vst.txt",
    targFile           = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10.txt",
    regFile            = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_final.txt",
    priorFile          = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv",
    priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"],
    taskLabelFile      = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/taskLabels.txt",  # NEW
    outputDir          = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/MultitaskInferelator/test",
    fusionLambda       = 0.1,           # fusion penalty (used by :fixed_ratio)
    similaritySource   = :expression,   # or :metadata / :ontology with similarityFile
    # ── Lambda optimization (pick one) ───────────────────────────────────────
    lambdaOpt          = :fixed_ratio,  # cheapest — good for first runs
    # lambdaOpt        = :ebic,         # principled — recommended for production
    # lambdaOpt        = :bstars_2d,    # strongest — use for benchmarking
    fusionRatio        = 0.1,           # for :fixed_ratio — lambda_f = fusionRatio * lambda_s
    ebicGamma          = 1.0,           # EBIC gamma — 1.0 recommended when p >> n
    gridSize           = 10,            # lambda grid size for :ebic and :bstars_2d
    refinementSize     = 10,            # fine grid size for :bstars_2d stage 2
    # ── Elastic net ──────────────────────────────────────────────────────────
    elasticNetAlpha    = 1.0            # 1.0=pure LASSO | 0.5=elastic net (recommended with demerged TFs)
)
