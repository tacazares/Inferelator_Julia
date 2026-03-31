# test/runtests.jl
using Test
using InferelatorJL
using Random
using Statistics
using LinearAlgebra

# =============================================================================
# Helpers — write synthetic input files to a temp directory
# =============================================================================

"""
Write a gene expression TSV:
  - Row 1: empty cell + sample names
  - Rows 2+: gene name + Float64 values
"""
function write_expr_tsv(path, geneNames, sampleNames, mat)
    open(path, "w") do io
        println(io, "\t" * join(sampleNames, "\t"))
        for (i, g) in enumerate(geneNames)
            println(io, g * "\t" * join(string.(mat[i, :]), "\t"))
        end
    end
end

"""
Write a gene list file (one name per line).
"""
write_gene_list(path, names) = write(path, join(names, "\n") * "\n")

"""
Write a sparse prior TSV:
  - Row 1: empty cell + TF names
  - Rows 2+: gene name + Float64 values
"""
function write_prior_tsv(path, tfNames, geneNames, mat)
    open(path, "w") do io
        println(io, "\t" * join(tfNames, "\t"))
        for (i, g) in enumerate(geneNames)
            println(io, g * "\t" * join(string.(mat[i, :]), "\t"))
        end
    end
end


# =============================================================================
@testset "InferelatorJL" begin
# =============================================================================


# -----------------------------------------------------------------------------
@testset "Structs instantiate" begin
# -----------------------------------------------------------------------------
    @test GeneExpressionData() isa GeneExpressionData
    @test PriorTFAData()       isa PriorTFAData
    @test mergedTFsResult()    isa mergedTFsResult
    @test GrnData()            isa GrnData
    @test BuildGrn()           isa BuildGrn
end


# -----------------------------------------------------------------------------
@testset "Data loading — expression matrix" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()

    # Genes intentionally out of alphabetical order to test sort
    genes   = ["Zap70", "Akt1", "Myc", "Foxp3"]
    samples = ["S1", "S2", "S3", "S4", "S5"]
    mat = Float64[
        1.0  2.0  3.0  4.0  5.0;   # Zap70
        6.0  7.0  8.0  9.0  10.0;  # Akt1
        0.1  0.2  0.1  0.2  0.1;   # Myc
        3.0  3.5  4.0  4.5  5.0;   # Foxp3
    ]
    exprFile = joinpath(tmpdir, "expr.txt")
    write_expr_tsv(exprFile, genes, samples, mat)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)

    # Gene names should be sorted alphabetically
    @test data.geneNames == sort(genes)

    # Matrix rows must match the sorted order
    sorted_order = sortperm(genes)        # [2,4,3,1] → Akt1,Foxp3,Myc,Zap70
    @test data.geneExpressionMat == mat[sorted_order, :]

    # Sample labels are parsed from the header
    @test data.cellLabels == samples

    # Dimensions: genes × samples
    @test size(data.geneExpressionMat) == (length(genes), length(samples))
end


# -----------------------------------------------------------------------------
@testset "Data loading — non-numeric first column values are gene names" begin
# -----------------------------------------------------------------------------
    # Verify that a gene name like "1500009L16Rik" (starts with a digit) is
    # treated as a gene name, not discarded or confused with a header.
    tmpdir  = mktempdir()
    genes   = ["1500009L16Rik", "Bcl2", "Cd44"]
    samples = ["A", "B", "C"]
    mat = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

    exprFile = joinpath(tmpdir, "expr_numnames.txt")
    write_expr_tsv(exprFile, genes, samples, mat)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)

    @test "1500009L16Rik" in data.geneNames
    @test length(data.geneNames) == 3
end


# -----------------------------------------------------------------------------
@testset "Data loading — invalid file path throws" begin
# -----------------------------------------------------------------------------
    data = GeneExpressionData()
    @test_throws ErrorException InferelatorJL.loadExpressionData!(data, "/nonexistent/path/expr.txt")
end


# -----------------------------------------------------------------------------
@testset "Data loading — target gene filtering" begin
# -----------------------------------------------------------------------------
    tmpdir  = mktempdir()
    genes   = ["Akt1", "Bcl2", "Myc", "Foxp3", "Stat3"]
    samples = ["S1", "S2", "S3", "S4", "S5", "S6"]

    # Myc has zero variance — should be removed with default epsilon
    mat = Float64[
        1.0  2.0  3.0  4.0  5.0  6.0;   # Akt1   — good variance
        2.0  2.5  3.0  2.5  2.0  3.5;   # Bcl2   — good variance
        1.0  1.0  1.0  1.0  1.0  1.0;   # Myc    — ZERO variance
        0.5  1.5  2.5  3.5  4.5  5.5;   # Foxp3  — good variance
        3.0  1.0  4.0  1.0  5.0  9.0;   # Stat3  — good variance
    ]
    exprFile = joinpath(tmpdir, "expr.txt")
    write_expr_tsv(exprFile, genes, samples, mat)

    # Request 4 genes (exclude Stat3); Myc will be removed by variance filter
    targFile = joinpath(tmpdir, "targs.txt")
    write_gene_list(targFile, ["Akt1", "Bcl2", "Myc", "Foxp3"])

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)

    # Myc must have been removed
    @test !("Myc" in data.targGenes)

    # The other three requested genes should be present
    @test "Akt1"  in data.targGenes
    @test "Bcl2"  in data.targGenes
    @test "Foxp3" in data.targGenes

    # targGeneMat dimensions must match retained genes × samples
    @test size(data.targGeneMat) == (3, length(samples))
end


# -----------------------------------------------------------------------------
@testset "Data loading — target file not found throws" begin
# -----------------------------------------------------------------------------
    tmpdir  = mktempdir()
    genes   = ["Akt1", "Bcl2"]
    samples = ["S1", "S2", "S3"]
    mat     = [1.0 2.0 3.0; 4.0 5.0 6.0]
    exprFile = joinpath(tmpdir, "expr.txt")
    write_expr_tsv(exprFile, genes, samples, mat)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(
        data, "/no/such/file.txt")
end


# -----------------------------------------------------------------------------
@testset "Data loading — potential regulators" begin
# -----------------------------------------------------------------------------
    tmpdir  = mktempdir()
    genes   = ["Akt1", "Bcl2", "Foxp3", "Myc"]
    samples = ["S1", "S2", "S3"]
    mat = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0; 10.0 11.0 12.0]
    exprFile = joinpath(tmpdir, "expr.txt")
    write_expr_tsv(exprFile, genes, samples, mat)

    # Regulator list: 2 present in expression, 1 not present
    regFile = joinpath(tmpdir, "regs.txt")
    write_gene_list(regFile, ["Foxp3", "Akt1", "NotInExpr"])

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    InferelatorJL.loadPotentialRegulators!(data, regFile)

    # Only the 2 regulators present in expression data should appear
    @test length(data.potRegs) == 2
    @test "Foxp3" in data.potRegs
    @test "Akt1"  in data.potRegs
    @test !("NotInExpr" in data.potRegs)

    # mRNA matrix: rows = regulators that have expression, cols = samples
    @test size(data.potRegMatmRNA) == (2, length(samples))
end


# =============================================================================
# Helper — build a fully loaded GeneExpressionData from scratch in memory
# =============================================================================
function make_test_data(tmpdir)
    genes   = ["Akt1", "Bcl2", "Foxp3", "Myc", "Stat3"]
    tfs     = ["Foxp3", "Stat3"]                # regulators
    samples = ["S$i" for i in 1:10]

    Random.seed!(42)
    mat = rand(Float64, length(genes), length(samples)) .* 5 .+ 0.5

    exprFile = joinpath(tmpdir, "expr.txt")
    targFile = joinpath(tmpdir, "targs.txt")
    regFile  = joinpath(tmpdir, "regs.txt")

    write_expr_tsv(exprFile, genes, samples, mat)
    write_gene_list(targFile, genes)           # all genes are targets
    write_gene_list(regFile,  tfs)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
    InferelatorJL.loadPotentialRegulators!(data, regFile)
    InferelatorJL.processTFAGenes!(data, "")   # use all genes for TFA

    return data, genes, tfs, samples
end


# -----------------------------------------------------------------------------
@testset "TFA — output dimensions" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()
    data, genes, tfs, samples = make_test_data(tmpdir)

    nTFs     = length(tfs)
    nGenes   = length(genes)
    nSamples = length(samples)

    # Prior: every gene has at least one edge per TF
    # Shape: genes × tfs (written as rows=genes, cols=tfs)
    prior_mat = Float64[
        1 0;   # Akt1  → Foxp3
        1 1;   # Bcl2  → Foxp3, Stat3
        0 1;   # Foxp3 → Stat3
        1 1;   # Myc   → both
        1 1;   # Stat3 → both
    ]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing,
                                    minTargets    = 1)
    InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

    # medTfas must be nTFs × nSamples  (result of priorMatrix \ targExpression)
    @test size(tfaData.medTfas, 1) == length(tfaData.pRegs)
    @test size(tfaData.medTfas, 2) == nSamples
end


# -----------------------------------------------------------------------------
@testset "TFA — least-squares correctness" begin
# -----------------------------------------------------------------------------
    # With a square, well-conditioned prior, TFA = prior \ expression exactly.
    tmpdir = mktempdir()

    genes   = ["Akt1", "Bcl2", "Foxp3"]
    tfs     = ["Foxp3", "Akt1"]
    samples = ["S1", "S2", "S3", "S4"]

    Random.seed!(7)
    expr_mat = rand(Float64, length(genes), length(samples)) .+ 1.0

    exprFile = joinpath(tmpdir, "expr.txt")
    targFile = joinpath(tmpdir, "targs.txt")
    regFile  = joinpath(tmpdir, "regs.txt")
    write_expr_tsv(exprFile, genes, samples, expr_mat)
    write_gene_list(targFile, genes)
    write_gene_list(regFile,  tfs)

    # Full prior: all 3 genes × 2 TFs, every cell nonzero (both TFs have >1 target)
    prior_mat = Float64[
        1.0  0.5;
        0.5  1.0;
        0.8  0.3;
    ]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, genes, prior_mat)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
    InferelatorJL.loadPotentialRegulators!(data, regFile)
    InferelatorJL.processTFAGenes!(data, "")

    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

    # Manually compute expected TFA: priorMatrix \ targExpression
    expected = tfaData.priorMatrix \ tfaData.targExpression
    @test tfaData.medTfas ≈ expected atol=1e-10
end


# -----------------------------------------------------------------------------
@testset "TFA — z-scored targets change output" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()
    data, genes, tfs, samples = make_test_data(tmpdir)

    prior_mat = Float64[
        1 0; 1 1; 0 1; 1 1; 1 1;
    ]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    tfaData_raw = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData_raw, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData_raw, data; edgeSS = 0, zTarget = false)

    tfaData_z = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData_z, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData_z, data; edgeSS = 0, zTarget = true)

    # z-scored and non-z-scored TFA should differ
    @test !isapprox(tfaData_raw.medTfas, tfaData_z.medTfas; atol=1e-6)
end


# -----------------------------------------------------------------------------
@testset "Prior — minTargets filter removes low-coverage TFs" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()
    data, genes, tfs, samples = make_test_data(tmpdir)

    # TF1 (Foxp3) has 4 targets; TF2 (Stat3) has only 1 target
    prior_mat = Float64[
        1 0;   # Akt1  → Foxp3 only
        1 0;   # Bcl2  → Foxp3 only
        1 1;   # Foxp3 → both
        1 0;   # Myc   → Foxp3 only
        0 1;   # Stat3 → Stat3 only (1 target)
    ]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    tfaData = PriorTFAData()
    # minTargets = 3: Foxp3 has 4 targets (passes), Stat3 has 1 (fails)
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 3)

    @test "Foxp3" in tfaData.pRegs
    @test !("Stat3" in tfaData.pRegs)
end


# -----------------------------------------------------------------------------
@testset "Penalty matrix — self-regulatory edges set to Inf in TFmRNA mode" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()
    data, genes, tfs, samples = make_test_data(tmpdir)

    prior_mat = Float64[
        1 0; 1 1; 0 1; 1 1; 1 1;
    ]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

    grnData = GrnData()
    # TFmRNA mode: self-regulatory edges must be Inf
    InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = "TFmRNA")
    InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                        priorFilePenalties = [priorFile],
                                        lambdaBias         = [0.5],
                                        tfaOpt             = "TFmRNA")

    # For each TF that is also a target gene, penalty[gene, TF] must be Inf
    for tf in data.potRegsmRNA
        targIdx = findfirst(==(tf), data.targGenes)
        tfIdx   = findfirst(==(tf), grnData.allPredictors)
        if targIdx !== nothing && tfIdx !== nothing
            @test grnData.penaltyMat[targIdx, tfIdx] == Inf
        end
    end
end


# -----------------------------------------------------------------------------
@testset "Subsamples — shape and reproducibility" begin
# -----------------------------------------------------------------------------
    tmpdir = mktempdir()
    data, genes, tfs, samples = make_test_data(tmpdir)

    prior_mat = Float64[1 0; 1 1; 0 1; 1 1; 1 1]
    priorFile = joinpath(tmpdir, "prior.tsv")
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

    grnData = GrnData()
    InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = "")

    totSS      = 10
    ssFrac     = 0.7
    nSamples   = length(data.cellLabels)
    expectCols = floor(Int, ssFrac * nSamples)

    Random.seed!(123)
    InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
    subsamps_a = copy(grnData.subsamps)

    # Correct shape
    @test size(subsamps_a) == (totSS, expectCols)

    # All indices are valid sample indices
    @test all(1 .<= subsamps_a .<= nSamples)

    # Reproducible with same seed
    Random.seed!(123)
    InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
    subsamps_b = copy(grnData.subsamps)
    @test subsamps_a == subsamps_b

    # Different seed → different subsamples (with overwhelming probability)
    Random.seed!(999)
    InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
    @test grnData.subsamps != subsamps_a
end


# -----------------------------------------------------------------------------
@testset "Edge cases — all target genes have zero variance" begin
# -----------------------------------------------------------------------------
    tmpdir  = mktempdir()
    genes   = ["Akt1", "Bcl2"]
    samples = ["S1", "S2", "S3"]
    # All rows constant — every gene will fail the variance filter
    mat = [1.0 1.0 1.0; 2.0 2.0 2.0]
    exprFile = joinpath(tmpdir, "expr.txt")
    targFile = joinpath(tmpdir, "targs.txt")
    write_expr_tsv(exprFile, genes, samples, mat)
    write_gene_list(targFile, genes)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(
        data, targFile; epsilon = 0.01)
end


# -----------------------------------------------------------------------------
@testset "Edge cases — target gene list requests genes not in expression data" begin
# -----------------------------------------------------------------------------
    tmpdir  = mktempdir()
    genes   = ["Akt1", "Bcl2"]
    samples = ["S1", "S2", "S3"]
    mat     = [1.0 2.0 3.0; 4.0 5.0 6.0]
    exprFile = joinpath(tmpdir, "expr.txt")
    targFile = joinpath(tmpdir, "targs.txt")
    write_expr_tsv(exprFile, genes, samples, mat)
    write_gene_list(targFile, ["NotAGene", "AlsoNotAGene"])

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(
        data, targFile; epsilon = 0.01)
end


# =============================================================================
end  # @testset "InferelatorJL"
# =============================================================================
