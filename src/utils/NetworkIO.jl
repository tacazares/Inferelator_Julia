using JLD2
using DelimitedFiles
using Printf
using Dates
using Statistics: std, median

# Save all core structs
function saveData(expressionData, tfaData, grnData, buildGrn, outputDir::String, fileName::String)
    output = joinpath(outputDir, fileName)
    @save output expressionData tfaData grnData buildGrn
end

"""
    writeBEBICLambdaTable!(grnData; outputDir)

Write a per-gene bEBIC lambda diagnostic TSV to `outputDir/bebic_lambda_summary.tsv`.

Columns:
- `gene`          : target gene name
- `median_lambda` : median EBIC-optimal lambda across subsamples
- `std_lambda`    : standard deviation of per-subsample EBIC-optimal lambdas —
                    high values indicate that EBIC is choosing very different
                    regularisation levels across subsamples for that gene

Requires `bebicEstimateLambdas!` to have been called first (populates `grnData.lambdaSS`).
"""
function writeBEBICLambdaTable!(grnData::GrnData; outputDir::String)
    genes         = grnData.targGenes
    medianLambdas = vec(median(grnData.lambdaSS, dims=2))
    stdLambdas    = [std(grnData.lambdaSS[i, :]) for i in 1:length(genes)]

    lines = ["gene\tmedian_lambda\tstd_lambda"]
    for i in 1:length(genes)
        push!(lines, "$(genes[i])\t$(medianLambdas[i])\t$(stdLambdas[i])")
    end

    outPath = joinpath(outputDir, "bebic_lambda_summary.tsv")
    open(outPath, "w") do io
        write(io, join(lines, "\n") * "\n")
    end
    @info "bEBIC lambda summary saved" path=outPath genes=length(genes)
end


"""
    writeEBICLambdaTable!(grnData; outputDir)

Write a per-gene EBIC lambda diagnostic TSV to `outputDir/ebic_lambda_summary.tsv`.

Columns:
- `gene`            : target gene name
- `chosen_lambda`   : EBIC-optimal lambda selected for that gene
- `n_nonzero_coefs` : number of non-zero LASSO coefficients at the chosen lambda —
                      i.e., the number of TFs retained as regulators of this gene

Requires `ebicSelect!` to have been called first (populates `grnData.ebicLambdas`
and `grnData.ebicNonzero`).
"""
function writeEBICLambdaTable!(grnData::GrnData; outputDir::String)
    genes         = grnData.targGenes
    chosenLambdas = grnData.ebicLambdas
    nNonzero      = grnData.ebicNonzero

    lines = ["gene\tchosen_lambda\tn_nonzero_coefs"]
    for i in 1:length(genes)
        push!(lines, "$(genes[i])\t$(chosenLambdas[i])\t$(nNonzero[i])")
    end

    outPath = joinpath(outputDir, "ebic_lambda_summary.tsv")
    open(outPath, "w") do io
        write(io, join(lines, "\n") * "\n")
    end
    @info "EBIC lambda summary saved" path=outPath genes=length(genes)
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