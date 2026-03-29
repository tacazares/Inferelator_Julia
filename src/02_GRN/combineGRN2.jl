"""
GRN

Main function:
- `combineGRNS2`: Combines GRNs using merged TFs, TFA data, and other metadata.  

Dependencies:
- `Data.GeneExpressionData` and `Data.PriorTFAData`
- `Prior.mergeDegenerateTFs`
"""
    # Access dependencies from the package; no need to include files again
    # using ..Data           # For GeneExpressionData
    # using ..PriorTFA       # For PriorTFAData
    # using ..MergeDegenerate

    # # Other standard packages
    # using PyPlot
    # using Statistics
    # using CSV
    # using DelimitedFiles
    # using JLD2
    # using NamedArrays

    # Export only the function users need
    # export combineGRNS2

function combineGRNS2(data::GeneExpressionData, mergedTFsData::mergedTFsResult, tfaGeneFile::Union{String, Nothing}, priorFile, edgeSS, minTargets,  
                        geneExprFile::Union{String, Nothing}=nothing, 
                        targetGeneFile::Union{String, Nothing}=nothing, 
                        potRegFile::Union{String, Nothing}=nothing;
                        outputDir::Union{String, Nothing}=nothing)
    if !isnothing(outputDir)
        mkpath(outputDir)
    end
                    
    # Ensure required expression data are available and non-empty
    requiredFields = [
        :tfaGenes,
        :tfaGeneMat,
        :potRegs,
        :potRegsmRNA,
        :potRegMatmRNA
    ]

    missingFields = [field for field in requiredFields if isempty(getfield(data, field))]

    if !isempty(missingFields)
        println("Missing or empty fields in data: ", missingFields)
        println("Generating required data by loading from input files...")
    
        requiredFiles = [
            ("Gene Expression File", geneExprFile),
            ("Target Gene File", targetGeneFile),
            ("Potential Regulators File", potRegFile),
        ]
    
        missingFiles = [name for (name, path) in requiredFiles if path === nothing || !isfile(path)]
        if !isempty(missingFiles)
            error("Cannot generate data. Missing input files: ", missingFiles)
        end
    
        data = GeneExpressionData()
        loadExpressionData!(data, geneExprFile)
        loadAndFilterTargetGenes!(data, targetGeneFile; epsilon=0.01)
        loadPotentialRegulators!(data, potRegFile)
        processTFAGenes!(data, tfaGeneFile; outputDir = outputDir)
    end


    # 2. Integrate prior information for TFA estimation
    tfaData = PriorTFAData()
    processPriorFile!(tfaData, data, priorFile; mergedTFsData, minTargets = minTargets);
    calculateTFA!(tfaData, data; edgeSS = edgeSS, outputDir = outputDir); 

    # Save TFA as a text file# Save median TFA if outputDir is specified
    if !isnothing(outputDir)
        namedMedTFA = NamedArray(tfaData.medTfas)
        setnames!(namedMedTFA, tfaData.pRegs, 1)
        setnames!(namedMedTFA, data.cellLabels, 2)

        outputfile = joinpath(outputDir, "TFA.txt")
        open(outputfile, "w") do io
            writedlm(io, permutedims(data.cellLabels))
            writedlm(io, namedMedTFA)
        end
    end

    return tfaData
end
