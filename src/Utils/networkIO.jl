using JLD2
using DelimitedFiles
using Printf
using Dates

# Save all core structs
function saveData(expressionData, tfaData, grnData, buildGrn, outputDir::String, fileName::String)
    output = joinpath(outputDir, fileName)
    @save output expressionData tfaData grnData buildGrn
end

# Write network tables
function writeNetworkTable!(buildGrn; outputDir::String, networkName::Union{String, Nothing}=nothing)
    baseName = (networkName === nothing || isempty(networkName)) ? "edges" : networkName * "edges"

    outputFile = joinpath(outputDir, baseName * ".tsv")
    colNames = "TF\tGene\tsignedQuantile\tStability\tCorrelation\tinPrior\n"
    open(outputFile, "w") do io
        write(io, colNames)
        writedlm(io, buildGrn.networkMat)
    end

    outputFileSubset = joinpath(outputDir, baseName * "_subset.tsv")
    open(outputFileSubset, "w") do io
        write(io, colNames)
        writedlm(io, buildGrn.networkMatSubset)
    end
end
