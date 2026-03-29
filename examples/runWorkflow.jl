cd("/data/miraldiNB/Michael/Scripts/GRN/InferelatorJL")
using Pkg
Pkg.activate(".")
using Revise
include("src/Inferelator.jl")
using .InferelatorJL

function runInferelator(;
    geneExprFile::String,
    targFile::String,
    regFile::String,
    priorFile::String,
    priorFilePenalties::Vector{String},
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
    instabilityLevel::String = "Network",
    useMeanEdgesPerGeneMode::Bool = true,
    combineOpt::String = "max",
    zTarget::Bool = true
)

    # --- Build output directory name
    subsamplePct = subsampleFrac * 100
    subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
    lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")
    networkBaseName = lowercase(instabilityLevel) * "Lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                      string(meanEdgesPerGene) * "tfsPerGene_" * "subsamplePCT" * subsampleStr
    dirOut = joinpath(outputDir, networkBaseName)
    mkpath(dirOut)

    println("=== Workflow Configuration ===")
    println("Output Directory: ", dirOut)
    println("GeneExpression File Used: ", geneExprFile)
    println("Processing prior file: ", priorFile)
    println("Prior files for penalties: ", priorFilePenalties)
    println("Lambda Bias used is: ", lambdaBias)
    println("Subsample Fraction: ", subsampleFrac)
    println("useMeanEdgesPerGeneMode: ", useMeanEdgesPerGeneMode)
    println("=============================")

    # ------- 1. Load expression data
    data = GeneExpressionData()
    loadExpressionData!(data, geneExprFile)
    loadAndFilterTargetGenes!(data, targFile; epsilon=0.01)
    loadPotentialRegulators!(data, regFile)
    processTFAGenes!(data, tfaGeneFile; outputDir=dirOut)

    # ------- 2. Merge degenerate TFs (optional)
    mergedTFsData = mergedTFsResult()
    mergeDegenerateTFs(mergedTFsData, priorFile; fileFormat=2)

    # ------- 3. Estimate TFA
    tfaData = PriorTFAData()
    processPriorFile!(tfaData, data, priorFile; mergedTFsData, minTargets=minTargets)
    calculateTFA!(tfaData, data; edgeSS=edgeSS, zscoreTargExp=zTarget, outputDir=dirOut)

    # ------- 4. Build GRN for each TFA option
    for tfaOpt in tfaOptions
        instabilitiesDir = tfaOpt == "" ? joinpath(dirOut, "TFA") : joinpath(dirOut, "TFmRNA")
        mkpath(instabilitiesDir)

        grnData = GrnData()
        preparePredictorMat!(grnData, data, tfaData, tfaOpt)
        preparePenaltyMatrix!(data, grnData, priorFilePenalties, lambdaBias, tfaOpt)
        constructSubsamples(data, grnData; totSS=bstarsTotSS, subsampleFrac=subsampleFrac)
        bstarsWarmStart(data, tfaData, grnData; minLambda=minLambda, maxLambda=maxLambda, 
                        totLambdasBstars=totLambdasBstars, targetInstability=targetInstability, zTarget=zTarget)
        constructSubsamples(data, grnData; totSS=totSS, subsampleFrac=subsampleFrac)
        bstartsEstimateInstability(grnData; totLambdas=totLambdas, instabilityLevel=instabilityLevel, 
                                   zTarget=zTarget, outputDir=instabilitiesDir)

        buildGrn = BuildGrn()
        chooseLambda!(grnData, buildGrn; instabilityLevel=instabilityLevel, targetInstability=targetInstability)
        rankEdges!(data, tfaData, grnData, buildGrn; useMeanEdgesPerGeneMode=useMeanEdgesPerGeneMode, 
                   meanEdgesPerGene=meanEdgesPerGene, correlationWeight=correlationWeight, outputDir=instabilitiesDir)
        writeNetworkTable!(buildGrn; outputDir=instabilitiesDir)
    end

    # ------- 5. Combine networks
    combinedNetDir = joinpath(dirOut, "Combined")
    nets2combine = [
        joinpath(dirOut, "TFA", "edges.tsv"),
        joinpath(dirOut, "TFmRNA", "edges.tsv")
    ]
    combineGRNs(nets2combine; combineOpt=combineOpt, meanEdgesPerGene=meanEdgesPerGene,
                useMeanEdgesPerGeneMode=useMeanEdgesPerGeneMode, saveDir=combinedNetDir, saveName="")

    # ------- 6. Re-estimate TFA for combined network
    netsCombinedSparse = joinpath(combinedNetDir, "combined_" * combineOpt * "_sp.tsv")
    combineGRNS2(data, mergedTFsData, tfaGeneFile, netsCombinedSparse, edgeSS, minTargets,
                 geneExprFile, targFile, regFile; outputDir=combinedNetDir)

end

# --- Run
runInferelator(
    geneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt",
    targFile     = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt",
    regFile      = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt",
    priorFile    = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv",
    priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"],
    outputDir    = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/test"
)