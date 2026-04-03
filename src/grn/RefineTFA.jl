"""
GRN

Main function:
- `refineTFA`: Combines GRNs using merged TFs, TFA data, and other metadata.  

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

function refineTFA(data::GeneExpressionData, mergedTFsData::mergedTFsResult;
                   priorFile::String                  = "",
                   tfaGeneFile::String                = "",
                   edgeSS::Int                        = 0,
                   minTargets::Int                    = 3,
                   zTarget::Bool                      = true,
                   timeLagFile::String                = "",
                   timeLag::Real                      = 0.0,
                   geneExprFile::String               = "",
                   targFile::String                   = "",
                   regFile::String                    = "",
                   outputDir::Union{String, Nothing}  = nothing)
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
        @info "Missing or empty fields in data — will reload from input files" missingFields
    
        requiredFiles = [
            ("Gene Expression File", geneExprFile),
            ("Target Gene File",     targFile),
            ("Potential Regulators File", regFile),
        ]

        missingFiles = [name for (name, path) in requiredFiles if isempty(path) || !isfile(path)]
        if !isempty(missingFiles)
            error("Cannot generate data. Missing input files: ", missingFiles)
        end

        data = GeneExpressionData()
        loadExpressionData!(data, geneExprFile)
        loadAndFilterTargetGenes!(data, targFile; epsilon=0.01)
        loadPotentialRegulators!(data, regFile)
        processTFAGenes!(data, tfaGeneFile; outputDir = outputDir)
    end


    # 2. Integrate prior information for TFA estimation
    tfaData = PriorTFAData()
    processPriorFile!(tfaData, data, priorFile; mergedTFsData, minTargets = minTargets)
    calculateTFA!(tfaData, data; edgeSS = edgeSS, zTarget = zTarget, outputDir = outputDir)
    if !isempty(timeLagFile)
        applyTimeLag!(tfaData, data, timeLagFile, timeLag)
    end

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
