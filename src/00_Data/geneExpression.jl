# 00_Data/geneExpression.jl
module Data

    export GeneExpressionData, loadExpressionData!, loadAndFilterTargetGenes!, loadPotentialRegulators!, processTFAGenes!

    using DelimitedFiles
    using Statistics
    using JLD2
    using CSV
    using Arrow
    using DataFrames

    # Initialize data structure

    mutable struct GeneExpressionData
        cellLabels::Vector{String}
        geneNames::Vector{String} # Full gene list
        geneExpressionMat::Matrix{Float64} # Full Expression Matrix
        potRegMatmRNA::Matrix{Float64}
        potRegs::Vector{String}
        potRegsmRNA::Vector{String}
        targGenes::Vector{String} # Target gene list
        targGeneMat::Matrix{Float64} # Filtered expression matrix
        tfaGenes::Vector{String}
        tfaGeneMat::Matrix{Float64}

        function GeneExpressionData()
            return new(
                [],
                [],
                Matrix{Float64}(undef, 0, 0),
                Matrix{Float64}(undef, 0, 0),
                [],
                [],
                [],
                Matrix{Float64}(undef, 0, 0),
                [],
                Matrix{Float64}(undef, 0, 0)
            )
        end
    end


    # Load Expression Data
    """
        loadExpressionData!(data::GeneExpressionData, geneExprFile)

    Loads gene expression data from a specified file and updates the `GeneExpressionData` object.

    # Arguments
    - `data::GeneExpressionData`: The data object to be populated with gene expression details.
    - `geneExprFile::String`: The path to the gene expression data file, which can be in `.arrow` format or tab-delimited text format.

    # Updates
    - `data.cellLabels`: A vector of sample condition names, derived from the column headers (excluding the first).
    - `data.targGenes`: A vector of gene names, derived from the first column.
    - `data.targGeneMat`: A matrix of expression values, converted to `Float64`.

    # Raises
    - An error if the gene expression file does not exist or if its path is invalid.

    # Notes
    - Assumes that the first column in the file contains gene names and subsequent columns contain cell labels and expression data.
    - For text files, data is sorted by gene names.
    """
    function loadExpressionData!(data::GeneExpressionData, geneExprFile)
        if isfile(geneExprFile)
            if endswith(geneExprFile, ".arrow")
                dfArrow = Arrow.Table(geneExprFile)
                df = deepcopy(DataFrame(dfArrow))
                dfArrow = nothing;
                cellLabels = names(df)[2:end]  # Assume first column is gene names
                geneNames = df[:, 1]  # Extract gene names
                ncounts = Matrix(df[:, 2:end])
                df = nothing
            else
                fid = open(geneExprFile);
                C = readdlm(fid, '\t', '\n');
                close(fid)
                cellLabels = C[1, :];
                cellLabels = filter(!isempty, cellLabels);
                C = C[2:end, :];
                inds = sortperm(C[:, 1]);
                C = C[inds, :];
                geneNames = C[:, 1];
                ncounts = C[:, 2:end]
            end
            ncounts = convert(Matrix{Float64}, ncounts);

            data.cellLabels = String.(cellLabels);
            data.geneNames = String.(geneNames);
            data.geneExpressionMat = ncounts;
            ncounts = nothing;
        else
            error("Expression data file path is invalid.")
        end
    end


    """
        loadAndFilterTargetGenes!(data::GeneExpressionData, targetGeneFile; eps=1E-10)

    Loads target genes from a file, filters them based on presence and variance in the existing gene expression data, and updates the `GeneExpressionData` object.

    # Arguments
    - `data::GeneExpressionData`: The data object to be updated with filtered target gene information.
    - `targetGeneFile::String`: The path to the file containing target gene names.
    - `eps::Float64`: A keyword argument specifying the variance threshold for filtering target genes. Default is `1E-10`.

    # Updates
    - `data.targGenes`: A filtered vector of target gene names that are present in the gene expression data and meet the variance threshold.
    - `data.targGeneMat`: A matrix of expression values for these filtered target genes.

    # Raises
    - An error if the target gene file does not exist.
    - An error if no target genes are found in the gene expression data.
    - An error if all target genes are filtered out due to low variance.

    # Notes
    - This function finds matching genes and filters them based on a variance threshold to ensure sufficient expression variability across samples.
    """
    function loadAndFilterTargetGenes!(data::GeneExpressionData, targetGeneFile; epsilon=1E-10)
        if isfile(targetGeneFile)
            # Load in target gene file
            fid = open(targetGeneFile)
            targetGenes = readdlm(fid, String)
            close(fid)
            
            # Find all geneNames that are in the targer gene file
            inds = findall(in(targetGenes), data.geneNames)
            if isempty(inds)
                error("No target genes found in expression data!")
            end
            
            # Filter target genes and expression matrix
            targGenes = data.geneNames[inds]
            targGeneMat = data.geneExpressionMat[inds, :]

            # Filter genes by minimum variance cutoff
            stds = std(targGeneMat, dims=2)
            keep = [index[1] for index in findall(stds .>= epsilon)]
            targGenesFilter = targGenes[keep]
            targGeneMatFilter = targGeneMat[keep, :]
            if isempty(keep)
                error("All target genes removed due to low variance!")
            else
                println(length(targGenesFilter), " target genes retained after filtering")
            end
            data.targGenes = targGenesFilter
            data.targGeneMat = targGeneMatFilter
        else
            error("Target gene file not found.")
        end
    end

    # Load Regulators
    """
        loadPotentialRegulators!(data::GeneExpressionData, potRegFile)

    Loads potential regulators and updates the `GeneExpressionData` object with their expression data.

    # Arguments
    - `data::GeneExpressionData`: The data object to be populated with potential regulator information.
    - `potRegFile::String`: The path to the file containing potential regulator names.

    # Updates
    - `data.potRegs`: List of all potential regulator names.
    - `data.potRegsmRNA`: List of regulator names with corresponding mRNA expression data.
    - `data.potRegMatmRNA`: Expression matrix for regulators with mRNA data.

    # Raises
    - An error if the potential regulators file does not exist.

    # Notes
    - Only the potential regulators present in the existing gene expression data are included.
    """
    function loadPotentialRegulators!(data::GeneExpressionData, potRegFile)
        if isfile(potRegFile)
            fid = open(potRegFile)
            potentialRegs = readdlm(fid, String)
            close(fid)
            potentialRegs = vec(potentialRegs)
            
            indsPotRegs = findall(in(data.geneNames), potentialRegs)
            potentialRegs = potentialRegs[indsPotRegs]

            inds = findall(in(potentialRegs), data.geneNames)
            potRegsmRNA = data.geneNames[inds]
            potRegMatmRNA = data.geneExpressionMat[inds, :]
            
            data.potRegs = potentialRegs
            data.potRegsmRNA = potRegsmRNA
            data.potRegMatmRNA = potRegMatmRNA
        else
            error("Potential regulators file not found.")
        end
    end



    # Process genes for TFA calculation
    """
        processTFAGenes(file, geneSc, nCounts)

    Processes TFA genes by retrieving their expression data from the expression matrix.

    # Arguments
    - `file::String`: The path to the TFA genes file. If the path is empty or invalid, all genes are used.
    - `geneNames::Vector{String}`: A vector of gene names from the expression data.
    - `geneExpressionMat::Matrix{Float64}`: A matrix containing expression values for the genes.

    # Returns
    - `tfaGenes::Vector{String}`: A vector of TFA gene names with expression data.
    - `tfaGeneMat::Matrix{Float64}`: A matrix of expression values for the TFA genes.

    # Notes
    - If the specified file does not exist, all genes are considered for TFA.
    """
    function processTFAGenes!(data::GeneExpressionData, tfaGeneFile::Union{String, Nothing}; outputDir::Union{String, Nothing}=nothing)
        if (tfaGeneFile !== nothing) && (tfaGeneFile != "") && isfile(tfaGeneFile)
            tfaGenes = readlines(tfaGeneFile)
        else
            tfaGenes = data.geneNames
        end

        inds = findall(in(tfaGenes), data.geneNames)
        data.tfaGenes = data.geneNames[inds]
        data.tfaGeneMat = data.geneExpressionMat[inds, :]

        if outputDir !== nothing && outputDir !== ""

            outputFile = joinpath(outputDir, "geneExprMat.jld")
            save_object(outputFile, data)

            # cellLabels = data.cellLabels
            # geneNames = data.geneNames
            # geneExpressionMat = data.geneExpressionMat
            # potRegs = data.potRegs
            # potRegsmRNA = data.potRegsmRNA
            # potRegMatmRNA = data.potRegMatmRNA
            # targGenes = data.targGenes
            # targGeneMat = data.targGeneMat
            # tfaGenes = data.tfaGenes
            # tfaGeneMat = data.tfaGeneMat
            # @save outputFile cellLabels geneNames geneExpressionMat potRegs potRegsmRNA potRegMatmRNA targGenes targGeneMat tfaGenes tfaGeneMat

        end
    end

end
