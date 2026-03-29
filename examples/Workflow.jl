# cd("/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Testing")
# include("../julia_fxns/InferelatorFunction1.jl")
# include("../julia_fxns/mergeDegeneratePriorTFs.jl")
# include("../julia_fxns/InferelatorFunction2.jl")
# # include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/InferelatorFunction3.jl")
# include("../julia_fxns/InferelatorFunction3.jl")
# include("../julia_fxns/InferelatorFunction4.jl")
# # include("../julia_fxns/InferelatorFunction5.jl")
# include("../julia_fxns/CombineTRN/combineGRN.jl")
# include("../julia_fxns/CombineTRN/combineGRN2.jl")

# Activate the folder where InferelatorJL/src is located
# cd("/data/miraldiNB/Michael/Scripts/GRN/InferelatorJL/src")
using Pkg
Pkg.activate("/data/miraldiNB/Michael/Scripts/GRN/InferelatorJL")
using Revise
include("src/Inferelator.jl")
Pkg.instantiate()
using .InferelatorJL
using .InferelatorJL: Data
using .InferelatorJL: MergeDegenerate
using .InferelatorJL: PriorTFA

outputDir = "/data/miraldiNB/Michael/projects/GRN/mCD4T_Wayman/Inferelator/test"
# outputDir = "/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/metaCells/SEACells"

tfaOptions = ["", "TFmRNA"]
totSS = 80
bstarsTotSS = 5
subsampleFrac = 0.68
minLambda = 0.01
maxLambda = 0.5
totLambdasBstars = 20
totLambdas = 40
targetInstability = 0.05
meanEdgesPerGene = 20
correlationWeight = 1
minTargets = 3
edgeSS = 0
lambdaBias = [0.5]
instabilityLevel = "Network" # options: Gene or Network
useMeanEdgesPerGeneMode = true
leaveOutSampleList = ""
combineOpt = "max"

# geneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt"
geneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_downsample_0.25M.txt"  # VST downSampled
# geneExprFile = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/geneExpression/SingleCell/Tfh10_scRNA_Cells_TMsRemoved_logNorm_Counts.arrow"
# geneExprFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/counts_final.txt"  # - spike control not ignored during batch correction
# geneExprFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/counts_VST_combat_final.txt"  # spike control ignored during batch correction
# geneExprFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/counts_VST_combat_0.15M_final.txt"
# geneExprFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/scNormCounts.arrow"
# geneExprFile = "/data/miraldiNB/Michael/hCD4T_Katko/SEAcells/metacell_log1p.txt"

# targFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Prior/SubsetPriors/all_targs.txt"
targFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt"

# regFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Prior/SubsetPriors/all_TFs.txt"
regFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt"

priorFile = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"
# priorFile = "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/baseGRN/Tfh10_base_GRN_dedup_normF.tsv"
# priorFile = "/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/FIMO_Cicero_b.tsv"
# priorFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/MotifScan5kbTSS_b.tsv"
# priorFile = "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/MaxATAC_5kbTSS_Trac_MicroC_combined.tsv"
# priorFile = "/data/miraldiNB/Michael/hCD4T_Katko/CellOracle/baseGRN/hCD4T_base_GRN_dedup_b.tsv"
# priorFile = "/Volumes/miraldiNB/Michael/hCD4T_Katko/dataBank/Priors/majorCellType/Naive/5kbTSS_FIMOp5_b.tsv"
priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/GRN/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"]
tfaGeneFile = ""

subsamplePct = subsampleFrac * 100
subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")
networkBaseName = lowercase(instabilityLevel)*"Lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                  string(meanEdgesPerGene) * "tfsPerGene_"  * "subsamplePCT" * subsampleStr
dirOut = joinpath(outputDir,networkBaseName)
mkpath(dirOut)

println("=== Workflow Configuration ===")
println("Output Directory:", dirOut)
println("GeneExpression File Used: ", geneExprFile)
println("Processing prior file: ", priorFile)
println("Prior files for penalties: ", priorFilePenalties)
println("lambda Bias used is:", lambdaBias)
println("Subsample Fraction: ", subsampleFrac)
println("useMeanEdgesPerGeneMode:", useMeanEdgesPerGeneMode )
println("=============================")

# ------- 1. Import gene expression data, list of regulators, list of target genes
data = GeneExpressionData()
# InferelatorJL.Data.GeneExpressionData()
# To see what fileds are available, use `fieldnames(typeof(data))`
loadExpressionData!(data, geneExprFile);
loadAndFilterTargetGenes!(data, targFile; epsilon=0.01);
loadPotentialRegulators!(data, regFile);
processTFAGenes!(data, tfaGeneFile ; outputDir = dirOut);

mergedTFsData = mergedTFsResult()
# This is an optional step. 
mergeDegenerateTFs(mergedTFsData, priorFile; fileFormat = 2);

#  --------2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) 
#  ---------- using prior-based TFA and TF mRNA levels
tfaData = PriorTFAData()
processPriorFile!(tfaData, data, priorFile; mergedTFsData, minTargets = minTargets);
calculateTFA!(tfaData, data; edgeSS = edgeSS, outputDir = dirOut);   # use this if saving output
# calculateTFA!(tfaData, data; edgeSS = edgeSS);   # use this if not saving output


for tfaOpt in tfaOptions 
    instabilitiesDir = tfaOpt == "" ? joinpath(dirOut, "TFA") : joinpath(dirOut, "TFmRNA")
    mkpath(instabilitiesDir)
    # ------- 3. estimate Instabilities
    grnData = GrnData()
    preparePredictorMat!(grnData, data, tfaData, tfaOpt);
    preparePenaltyMatrix!(data, grnData, priorFilePenalties, lambdaBias, tfaOpt);
    constructSubsamples(data, grnData; totSS = bstarsTotSS, subsampleFrac = subsampleFrac);
    bstarsWarmStart(data, tfaData, grnData; minLambda = minLambda, maxLambda = maxLambda, totLambdasBstars = totLambdasBstars, targetInstability = targetInstability);
    constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = subsampleFrac);
    bstartsEstimateInstability(grnData; totLambdas = totLambdas, instabilityLevel = instabilityLevel, outputDir =  instabilitiesDir)

    # -------- 4. For a given instability cutoff and model size, rank TF-geneinteractions, and calculate stabilities 
    grnDir = instabilitiesDir
    buildGrn = BuildGrn()
    chooseLambda!(grnData, buildGrn; instabilityLevel = instabilityLevel, targetInstability = targetInstability)
    rankEdges!(data, tfaData, grnData, buildGrn; useMeanEdgesPerGeneMode = useMeanEdgesPerGeneMode, meanEdgesPerGene = meanEdgesPerGene, correlationWeight = correlationWeight, outputDir = grnDir)
    # saveData(data, tfaData, grnData, buildGrn,"/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/NewPipelineComparisonSc", "test_" * tfaOpt * ".jld")
    writeNetworkTable!(buildGrn; outputDir = grnDir)
end

# ------- 5. Combine networks for current prior and lambda.
combinedNetDir = joinpath(dirOut, "Combined")
nets2combine   = [
    joinpath(dirOut, "TFA", "edges.tsv"),
    joinpath(dirOut, "TFmRNA", "edges.tsv")
]
# combineGRNs(meanEdgesPerGene, useMeanEdgesPerGeneMode, combineOpt, combinedNetDir, nets2combine)
combineGRNs(nets2combine; combineOpt=combineOpt,
        meanEdgesPerGene=meanEdgesPerGene,
        useMeanEdgesPerGeneMode=useMeanEdgesPerGeneMode,
        saveDir=combinedNetDir,
        saveName=""
        )
netsCombinedSparse = joinpath(combinedNetDir,"combined_" * combineOpt* "_sp.tsv")

combineGRNS2(data, mergedTFsData, tfaGeneFile, netsCombinedSparse, edgeSS, minTargets,
                geneExprFile,targFile, regFile;
                outputDir=combinedNetDir)


# # # /Combined/combined_max.tsv"
# # mATACTracMicroTSS
# old =CSV.read("/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/Bulk/5kbTSS/newPipe/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.tsv", DataFrame; delim = "\t")
# new =CSV.read("/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/Bulk/5kbTSS/newPipe/spikeIgnored/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.tsv", DataFrame; delim = "\t")
# # # old2 = CSV.read("/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/5kbTSS/Bulk/oldPrior/geneLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", DataFrame; delim = "\t")
# # # new2 = CSV.read("/data/miraldiNB/Michael/hCD4T_Katko/Inferelator/noMergedTF/5kbTSS/Bulk/networkLambda0p5_100totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", DataFrame; delim = "\t")

# newPairs = Set((row.TF, row.Gene) for row in eachrow(new))
# oldPairs = Set((row.TF, row.Gene) for row in eachrow(old))
# # # old2Pairs = Set((row.TF, row.Gene) for row in eachrow(old2))
# # # new2Pairs = Set((row.TF, row.Gene) for row in eachrow(new2))

# # # intersectPairs = length(intersect(oldPairs, old2Pairs))
# intersectPairs = length(intersect(newPairs, oldPairs))
# # # intersectPairs = length(intersect(newPairs, old2Pairs))
# # # intersectPairs = length(intersect(new2Pairs, oldPairs))
# # # intersectPairs = length(intersect(new2Pairs, old2Pairs))
# # # length(intersect(new2Pairs, newPairs))






# ct <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/counts_VST_combat_final.txt", header = T)
# # Naive
# Naive <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/Naive.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # TCM
# TCM <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/TCM.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# colnames(TCM) <- gsub("\\.", "_", colnames(TCM))
# # TEM
# TEM <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/TEM.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# colnames(TEM) <- gsub("\\.", "_", colnames(TEM))
# # Treg
# Treg <- read.table("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/Treg.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# colnames(Treg) <- gsub("\\.", "_", colnames(Treg))


# dirOut <- "/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/spikeIgnored"
# dir.create("/data/miraldiNB/Michael/hCD4T_Katko/dataBank/GEX/majorCellType/spikeIgnored", recursive = T, showWarnings =T )

# Naive1 <- ct[ , setdiff(colnames(Naive), "X")]
# TCM1 <- ct[ , setdiff(colnames(TCM), "X")]
# TEM1 <- ct[ , setdiff(colnames(TEM), "X")]
# Treg1 <- ct[ , setdiff(colnames(Treg), "X")]

# write.table(Treg1, file.path(dirOut, "TREG.tsv"), col.names = NA, row.names = T, sep = "\t", quote = F)
