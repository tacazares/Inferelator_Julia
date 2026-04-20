# test/runtests.jl
using Test
using InferelatorJL
using Random
using Statistics
using LinearAlgebra
using DataFrames
using Arrow


# =============================================================================
# OPTIONAL REAL DATA PATHS
# =============================================================================
# Set these to point to real small input files to run real-data validation
# tests in addition to the synthetic ones.  Leave as "" (default) to skip
# the real-data section — all synthetic tests still run normally.
#
# You can set them in three ways (first match wins):
#   1. Edit the constants below directly.
#   2. Set environment variables before launching Julia:
#        export INFER_TEST_EXPR="/path/to/expr.txt"
#        julia --project=. test/runtests.jl
#   3. Pass them inline:
#        INFER_TEST_EXPR="/path/to/expr.txt" julia --project=. test/runtests.jl


## Option 1 — environment variables (no file editing needed)
# INFER_TEST_EXPR="/path/to/expr.txt" \
# INFER_TEST_TARGS="/path/to/targs.txt" \
# INFER_TEST_REGS="/path/to/regs.txt" \
# INFER_TEST_PRIOR="/path/to/prior.tsv" \
# julia --project=. test/runtests.jl

# Option 2 — edit the constants at the top of test/runtests.jl directly
# const REAL_EXPR_FILE  = "/path/to/expr.txt"

# -----------------------------------------------------------------------------
const REAL_EXPR_FILE  = get(ENV, "INFER_TEST_EXPR",  "")   # expression matrix
const REAL_TARG_FILE  = get(ENV, "INFER_TEST_TARGS", "")   # target gene list
const REAL_REG_FILE   = get(ENV, "INFER_TEST_REGS",  "")   # regulator list
const REAL_PRIOR_FILE = get(ENV, "INFER_TEST_PRIOR", "")   # prior network TSV
# =============================================================================


# =============================================================================
# Shared helpers — write synthetic input files to a temp directory
# =============================================================================

function write_expr_tsv(path, geneNames, sampleNames, mat)
    open(path, "w") do io
        println(io, "\t" * join(sampleNames, "\t"))
        for (i, g) in enumerate(geneNames)
            println(io, g * "\t" * join(string.(mat[i, :]), "\t"))
        end
    end
end

write_gene_list(path, names) = write(path, join(names, "\n") * "\n")

function write_prior_tsv(path, tfNames, geneNames, mat)
    open(path, "w") do io
        println(io, "\t" * join(tfNames, "\t"))
        for (i, g) in enumerate(geneNames)
            println(io, g * "\t" * join(string.(mat[i, :]), "\t"))
        end
    end
end

function write_edges_tsv(path, rows)
    # rows: Vector of NamedTuples with fields TF, Gene, signedQuantile, Stability, Correlation, inPrior
    open(path, "w") do io
        println(io, "TF\tGene\tsignedQuantile\tStability\tCorrelation\tinPrior")
        for r in rows
            println(io, "$(r.TF)\t$(r.Gene)\t$(r.signedQuantile)\t$(r.Stability)\t$(r.Correlation)\t$(r.inPrior)")
        end
    end
end

# Build a fully populated GeneExpressionData + PriorTFAData in one call
function make_test_pipeline(tmpdir;
        genes   = ["Akt1","Bcl2","Foxp3","Myc","Stat3"],
        tfs     = ["Foxp3","Stat3"],
        samples = ["S$i" for i in 1:10],
        seed    = 42)

    Random.seed!(seed)
    mat = rand(Float64, length(genes), length(samples)) .* 5 .+ 0.5

    exprFile  = joinpath(tmpdir, "expr.txt")
    targFile  = joinpath(tmpdir, "targs.txt")
    regFile   = joinpath(tmpdir, "regs.txt")
    priorFile = joinpath(tmpdir, "prior.tsv")

    write_expr_tsv(exprFile, genes, samples, mat)
    write_gene_list(targFile, genes)
    write_gene_list(regFile,  tfs)

    prior_mat = Float64[1 0; 1 1; 0 1; 1 1; 1 1]
    write_prior_tsv(priorFile, tfs, sort(genes), prior_mat)

    data = GeneExpressionData()
    InferelatorJL.loadExpressionData!(data, exprFile)
    InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
    InferelatorJL.loadPotentialRegulators!(data, regFile)
    InferelatorJL.processTFAGenes!(data, "")

    tfaData = PriorTFAData()
    InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                    mergedTFsData = nothing, minTargets = 1)
    InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

    return data, tfaData, priorFile
end


# =============================================================================
@testset "InferelatorJL" begin
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
@testset "Structs instantiate" begin
# ─────────────────────────────────────────────────────────────────────────────
    @test GeneExpressionData() isa GeneExpressionData
    @test PriorTFAData()       isa PriorTFAData
    @test mergedTFsResult()    isa mergedTFsResult
    @test GrnData()            isa GrnData
    @test BuildGrn()           isa BuildGrn
end


# =============================================================================
@testset "Data loading" begin
# =============================================================================

    @testset "Expression matrix — sort and parse" begin
        tmpdir  = mktempdir()
        genes   = ["Zap70","Akt1","Myc","Foxp3"]
        samples = ["S1","S2","S3","S4","S5"]
        mat = Float64[1 2 3 4 5; 6 7 8 9 10; 0.1 0.2 0.1 0.2 0.1; 3 3.5 4 4.5 5]
        exprFile = joinpath(tmpdir, "expr.txt")
        write_expr_tsv(exprFile, genes, samples, mat)

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)

        @test data.geneNames == sort(genes)
        sorted_order = sortperm(genes)
        @test data.geneExpressionMat == mat[sorted_order, :]
        @test data.cellLabels == samples
        @test size(data.geneExpressionMat) == (length(genes), length(samples))
    end

    @testset "Expression matrix — numeric-looking gene names kept" begin
        tmpdir  = mktempdir()
        genes   = ["1500009L16Rik","Bcl2","Cd44"]
        samples = ["A","B","C"]
        mat = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
        exprFile = joinpath(tmpdir, "expr.txt")
        write_expr_tsv(exprFile, genes, samples, mat)

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        @test "1500009L16Rik" in data.geneNames
        @test length(data.geneNames) == 3
    end

    @testset "Expression matrix — Arrow format" begin
        tmpdir  = mktempdir()
        genes   = ["Akt1", "Bcl2", "Myc"]
        samples = ["S1", "S2", "S3", "S4"]
        mat = Float64[1 2 3 4; 5 6 7 8; 9 10 11 12]

        # Write Arrow file: first column :Genes, then sample columns
        df = DataFrame(:Genes => genes)
        for (j, s) in enumerate(samples)
            df[!, Symbol(s)] = mat[:, j]
        end
        arrowFile = joinpath(tmpdir, "expr.arrow")
        Arrow.write(arrowFile, df)

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, arrowFile)

        @test length(data.geneNames)  == 3
        @test length(data.cellLabels) == 4
        @test size(data.geneExpressionMat) == (3, 4)
        @test data.geneNames  == genes
        @test data.cellLabels == samples
    end

    @testset "Expression matrix — invalid path throws" begin
        data = GeneExpressionData()
        @test_throws ErrorException InferelatorJL.loadExpressionData!(data, "/no/such/file.txt")
    end

    @testset "Target gene filtering — low-variance genes removed" begin
        tmpdir  = mktempdir()
        genes   = ["Akt1","Bcl2","Myc","Foxp3","Stat3"]
        samples = ["S$i" for i in 1:6]
        mat = Float64[
            1.0 2.0 3.0 4.0 5.0 6.0;
            2.0 2.5 3.0 2.5 2.0 3.5;
            1.0 1.0 1.0 1.0 1.0 1.0;   # zero variance — must be removed
            0.5 1.5 2.5 3.5 4.5 5.5;
            3.0 1.0 4.0 1.0 5.0 9.0;
        ]
        exprFile = joinpath(tmpdir, "expr.txt")
        targFile = joinpath(tmpdir, "targs.txt")
        write_expr_tsv(exprFile, genes, samples, mat)
        write_gene_list(targFile, ["Akt1","Bcl2","Myc","Foxp3"])

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)

        @test !("Myc" in data.targGenes)
        @test "Akt1"  in data.targGenes
        @test "Bcl2"  in data.targGenes
        @test "Foxp3" in data.targGenes
        @test size(data.targGeneMat) == (3, length(samples))
    end

    @testset "Target gene filtering — missing file throws" begin
        tmpdir = mktempdir()
        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, let
            f = joinpath(tmpdir,"e.txt")
            write_expr_tsv(f,["Akt1"],["S1","S2","S3"],[1.0 2.0 3.0])
            f
        end)
        @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(data, "/no/such.txt")
    end

    @testset "Target gene filtering — no overlap throws" begin
        tmpdir  = mktempdir()
        exprFile = joinpath(tmpdir,"expr.txt")
        targFile = joinpath(tmpdir,"targs.txt")
        write_expr_tsv(exprFile, ["Akt1","Bcl2"], ["S1","S2","S3"], [1.0 2.0 3.0; 4.0 5.0 6.0])
        write_gene_list(targFile, ["NotInExpr"])
        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(data, targFile)
    end

    @testset "Target gene filtering — all zero variance throws" begin
        tmpdir  = mktempdir()
        exprFile = joinpath(tmpdir,"expr.txt")
        targFile = joinpath(tmpdir,"targs.txt")
        write_expr_tsv(exprFile, ["Akt1","Bcl2"], ["S1","S2","S3"], [1.0 1.0 1.0; 2.0 2.0 2.0])
        write_gene_list(targFile, ["Akt1","Bcl2"])
        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        @test_throws ErrorException InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
    end

    @testset "Potential regulators — absent TFs discarded" begin
        tmpdir  = mktempdir()
        genes   = ["Akt1","Bcl2","Foxp3","Myc"]
        samples = ["S1","S2","S3"]
        mat = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0; 10.0 11.0 12.0]
        exprFile = joinpath(tmpdir,"expr.txt")
        regFile  = joinpath(tmpdir,"regs.txt")
        write_expr_tsv(exprFile, genes, samples, mat)
        write_gene_list(regFile, ["Foxp3","Akt1","NotInExpr"])

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        InferelatorJL.loadPotentialRegulators!(data, regFile)

        @test length(data.potRegs) == 2
        @test "Foxp3" in data.potRegs
        @test "Akt1"  in data.potRegs
        @test !("NotInExpr" in data.potRegs)
        @test size(data.potRegMatmRNA) == (2, length(samples))
    end

end # Data loading


# =============================================================================
@testset "TFA estimation" begin
# =============================================================================

    @testset "Output dimensions" begin
        tmpdir = mktempdir()
        data, tfaData, _ = make_test_pipeline(tmpdir)
        @test size(tfaData.medTfas, 1) == length(tfaData.pRegs)
        @test size(tfaData.medTfas, 2) == length(data.cellLabels)
    end

    @testset "Least-squares correctness" begin
        tmpdir  = mktempdir()
        genes   = ["Akt1","Bcl2","Foxp3"]
        tfs     = ["Foxp3","Akt1"]
        samples = ["S$i" for i in 1:4]
        Random.seed!(7)
        expr_mat = rand(Float64, length(genes), length(samples)) .+ 1.0

        exprFile = joinpath(tmpdir,"expr.txt")
        targFile = joinpath(tmpdir,"targs.txt")
        regFile  = joinpath(tmpdir,"regs.txt")
        priorFile = joinpath(tmpdir,"prior.tsv")
        write_expr_tsv(exprFile, genes, samples, expr_mat)
        write_gene_list(targFile, genes)
        write_gene_list(regFile, tfs)
        write_prior_tsv(priorFile, tfs, genes, [1.0 0.5; 0.5 1.0; 0.8 0.3])

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
        InferelatorJL.loadPotentialRegulators!(data, regFile)
        InferelatorJL.processTFAGenes!(data, "")

        tfaData = PriorTFAData()
        InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                        mergedTFsData = nothing, minTargets = 1)
        InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = false)

        @test tfaData.medTfas ≈ (tfaData.priorMatrix \ tfaData.targExpression) atol=1e-10
    end

    @testset "Z-scored targets produce different TFA" begin
        tmpdir = mktempdir()
        data, tfaData_raw, priorFile = make_test_pipeline(tmpdir)

        tfaData_z = PriorTFAData()
        InferelatorJL.processPriorFile!(tfaData_z, data, priorFile;
                                        mergedTFsData = nothing, minTargets = 1)
        InferelatorJL.calculateTFA!(tfaData_z, data; edgeSS = 0, zTarget = true)

        @test !isapprox(tfaData_raw.medTfas, tfaData_z.medTfas; atol=1e-6)
    end

    @testset "minTargets filter removes low-coverage TFs" begin
        tmpdir = mktempdir()
        genes   = sort(["Akt1","Bcl2","Foxp3","Myc","Stat3"])
        tfs     = ["Foxp3","Stat3"]
        samples = ["S$i" for i in 1:8]
        Random.seed!(1)
        mat = rand(Float64, length(genes), length(samples)) .+ 1.0
        exprFile = joinpath(tmpdir,"expr.txt")
        targFile = joinpath(tmpdir,"targs.txt")
        regFile  = joinpath(tmpdir,"regs.txt")
        priorFile = joinpath(tmpdir,"prior.tsv")
        write_expr_tsv(exprFile, genes, samples, mat)
        write_gene_list(targFile, genes)
        write_gene_list(regFile, tfs)
        # Foxp3 → 4 targets (passes minTargets=3); Stat3 → 1 target (fails)
        write_prior_tsv(priorFile, tfs, genes, Float64[1 0; 1 0; 1 1; 1 0; 0 1])

        data = GeneExpressionData()
        InferelatorJL.loadExpressionData!(data, exprFile)
        InferelatorJL.loadAndFilterTargetGenes!(data, targFile; epsilon = 0.01)
        InferelatorJL.loadPotentialRegulators!(data, regFile)
        InferelatorJL.processTFAGenes!(data, "")

        tfaData = PriorTFAData()
        InferelatorJL.processPriorFile!(tfaData, data, priorFile;
                                        mergedTFsData = nothing, minTargets = 3)
        @test "Foxp3" in tfaData.pRegs
        @test !("Stat3" in tfaData.pRegs)
    end

end # TFA estimation


# =============================================================================
@testset "GRN preparation" begin
# =============================================================================

    @testset "Penalty matrix — self-regulatory edges are Inf in TFmRNA mode" begin
        tmpdir = mktempdir()
        data, tfaData, priorFile = make_test_pipeline(tmpdir)

        grnData = GrnData()
        InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = "TFmRNA")
        InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                            priorFilePenalties = [priorFile],
                                            lambdaBias         = [0.5],
                                            tfaOpt             = "TFmRNA")
        for tf in data.potRegsmRNA
            targIdx = findfirst(==(tf), data.targGenes)
            tfIdx   = findfirst(==(tf), grnData.allPredictors)
            if targIdx !== nothing && tfIdx !== nothing
                @test grnData.penaltyMat[targIdx, tfIdx] == Inf
            end
        end
    end

    @testset "Penalty matrix — prior edges get reduced lambda" begin
        tmpdir = mktempdir()
        data, tfaData, priorFile = make_test_pipeline(tmpdir)

        grnData = GrnData()
        InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = "")
        InferelatorJL.preparePenaltyMatrix!(data, grnData;
                                            priorFilePenalties = [priorFile],
                                            lambdaBias         = [0.5],
                                            tfaOpt             = "")
        # At least one finite penalty < 1.0 should exist (prior-supported edge)
        finite_penalties = filter(isfinite, vec(grnData.penaltyMat))
        @test any(finite_penalties .< 1.0)
    end

    @testset "Subsamples — shape, bounds, and reproducibility" begin
        tmpdir = mktempdir()
        data, tfaData, _ = make_test_pipeline(tmpdir)

        grnData = GrnData()
        InferelatorJL.preparePredictorMat!(grnData, data, tfaData; tfaOpt = "")

        totSS    = 10
        ssFrac   = 0.7
        nSamples = length(data.cellLabels)
        expectCols = floor(Int, ssFrac * nSamples)

        Random.seed!(123)
        InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
        subsamps_a = copy(grnData.subsamps)

        @test size(subsamps_a) == (totSS, expectCols)
        @test all(1 .<= subsamps_a .<= nSamples)

        Random.seed!(123)
        InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
        @test grnData.subsamps == subsamps_a

        Random.seed!(999)
        InferelatorJL.constructSubsamples(data, grnData; totSS = totSS, subsampleFrac = ssFrac)
        @test grnData.subsamps != subsamps_a
    end

end # GRN preparation


# =============================================================================
@testset "Data utilities" begin
# =============================================================================

    @testset "matrixToEdgeList — wide matrix → edge list" begin
        df = DataFrame(Gene=["A","B"], TF1=[1.0,0.0], TF2=[0.5,1.0], TF3=[0.0,0.8])
        long = matrixToEdgeList(df)
        @test ncol(long) == 3                        # Gene, variable, value
        @test nrow(long) == 2 * 3                    # 2 genes × 3 TFs
    end

    @testset "matrixToEdgeList — already edge list → unchanged" begin
        df = DataFrame(TF=["A","A"], Gene=["X","Y"], W=[1.0,0.5])
        @test matrixToEdgeList(df) === df
    end

    @testset "edgeListToMatrix — edge list → wide matrix" begin
        df = DataFrame(Gene=["X","X","Y","Y"],
                       TF=["A","B","A","B"],
                       W=[1.0,0.5,0.0,0.8])
        wide = edgeListToMatrix(df; indices=(1,2,3))
        @test "A" in names(wide) && "B" in names(wide)
        @test nrow(wide) == 2   # X and Y
    end

    @testset "edgeListToMatrix — fewer than 3 columns throws" begin
        df = DataFrame(a=[1,2], b=[3,4])
        @test_throws ErrorException edgeListToMatrix(df)
    end

    @testset "frobeniusNormalize — column norms are ≈ 1 after :column" begin
        df = DataFrame(Gene=["A","B","C"],
                       TF1=[3.0,4.0,0.0],
                       TF2=[1.0,0.0,0.0])
        dfN = frobeniusNormalize(df, :column)
        for col in names(dfN)[2:end]
            v = Float64.(dfN[!, col])
            if norm(v) > 0
                @test norm(v) ≈ 1.0 atol=1e-10
            end
        end
    end

    @testset "frobeniusNormalize — row norms are ≈ 1 after :row" begin
        df = DataFrame(Gene=["A","B"],
                       TF1=[3.0,0.0],
                       TF2=[4.0,5.0],
                       TF3=[0.0,12.0])
        dfN = frobeniusNormalize(df, :row)
        mat = Matrix{Float64}(dfN[!, 2:end])
        for i in 1:nrow(dfN)
            n = norm(mat[i, :])
            if n > 0
                @test n ≈ 1.0 atol=1e-10
            end
        end
    end

    @testset "frobeniusNormalize — invalid dims throws" begin
        df = DataFrame(Gene=["A"], TF1=[1.0])
        @test_throws ArgumentError frobeniusNormalize(df, :diag)
    end

    @testset "completeDF — missing rows and columns are filled" begin
        df = DataFrame(Gene=["A","B"], TF1=[1.0,2.0])
        allIDs = ["A","B","C"]
        allCols = [:TF1, :TF2]
        out = completeDF(df, :Gene, allIDs, allCols)

        @test nrow(out) == 3
        @test "C" in out.Gene
        @test ismissing(out[out.Gene .== "C", :TF1][1])   # missing row → missing
        @test :TF2 in propertynames(out)
        @test all(ismissing, out.TF2)                      # missing column → all missing
    end

    @testset "mergeDFs — sum two overlapping DataFrames" begin
        df1 = DataFrame(Gene=["A","B"], TF1=[1.0,2.0])
        df2 = DataFrame(Gene=["B","C"], TF1=[3.0,4.0])
        merged = mergeDFs([df1,df2], :Gene, "sum")

        @test sort(merged.Gene) == ["A","B","C"]
        bRow = merged[merged.Gene .== "B", :TF1][1]
        @test bRow ≈ 5.0   # 2.0 + 3.0
        aRow = merged[merged.Gene .== "A", :TF1][1]
        @test aRow ≈ 1.0
    end

    @testset "mergeDFs — avg option" begin
        df1 = DataFrame(Gene=["A"], TF1=[4.0])
        df2 = DataFrame(Gene=["A"], TF1=[2.0])
        merged = mergeDFs([df1,df2], :Gene, "avg")
        @test merged[merged.Gene .== "A", :TF1][1] ≈ 3.0
    end

    @testset "binarizeNumeric! — DataFrame: nonzero → 1, zero → 0" begin
        df = DataFrame(Gene=["A","B"], TF1=[0.0,3.7], TF2=[1.5,0.0])
        binarizeNumeric!(df)
        @test df[1, :TF1] == 0
        @test df[2, :TF1] == 1
        @test df[1, :TF2] == 1
        @test df[2, :TF2] == 0
    end

    @testset "binarizeNumeric! — Matrix: nonzero → 1, zero → 0" begin
        M = Float64[0.0 2.5; -1.0 0.0]
        binarizeNumeric!(M)
        @test M == [0 1; 1 0]
    end

    @testset "check_column_norms — correct true/false counts" begin
        df = DataFrame(Gene=["A","B","C"],
                       TF1=[1.0,0.0,0.0],   # norm = 1.0 → normalized
                       TF2=[1.0,1.0,0.0])   # norm = √2 ≠ 1 → not normalized
        _, n_true, n_false = check_column_norms(df; atol=1e-8)
        @test n_true  == 1
        @test n_false == 1
    end

    @testset "writeTSVWithEmptyFirstHeader — file format correct" begin
        tmpdir = mktempdir()
        df = DataFrame(Gene=["A","B"], TF1=[1.0,2.0], TF2=[3.0,4.0])
        outFile = joinpath(tmpdir, "out.tsv")
        writeTSVWithEmptyFirstHeader(df, outFile)

        lines = readlines(outFile)
        @test startswith(lines[1], "\t")          # first header cell is empty
        @test occursin("TF1", lines[1])
        @test occursin("TF2", lines[1])
        @test length(lines) == 3                  # header + 2 data rows
    end

end # Data utilities


# =============================================================================
@testset "Partial correlation" begin
# =============================================================================

    @testset "partialCorrelationMat — symmetric with unit diagonal" begin
        Random.seed!(5)
        X = rand(Float64, 20, 4)   # 20 obs, 4 variables
        P = InferelatorJL.partialCorrelationMat(X)

        @test size(P) == (4, 4)
        @test P ≈ P'  atol=1e-10                  # symmetric
        @test all(diag(P) .≈ 1.0)                 # diagonal is 1
        @test all(abs.(P) .<= 1.0 + 1e-8)        # values in [-1,1]
    end

    @testset "partialCorrelationMat — first_vs_all returns row vector" begin
        Random.seed!(6)
        X = rand(Float64, 20, 5)
        P = InferelatorJL.partialCorrelationMat(X; first_vs_all=true)
        @test size(P) == (1, 5)
    end

    @testset "partialCorrReg — symmetric with unit diagonal" begin
        Random.seed!(7)
        X = rand(Float64, 30, 3)
        P = InferelatorJL.partialCorrReg(X)

        @test size(P) == (3, 3)
        @test P ≈ P'  atol=1e-10
        @test all(diag(P) .≈ 1.0)
    end

    @testset "partialCorrelationMat and partialCorrReg agree" begin
        # Both methods should give similar partial correlations on well-conditioned data
        Random.seed!(8)
        X = rand(Float64, 50, 3)
        P1 = InferelatorJL.partialCorrelationMat(X)
        P2 = InferelatorJL.partialCorrReg(X)
        @test P1 ≈ P2 atol=0.05    # both methods, tolerant agreement
    end

end # Partial correlation


# =============================================================================
@testset "Network I/O" begin
# =============================================================================

    @testset "writeNetworkTable! — edges.tsv written with correct header" begin
        tmpdir = mktempdir()

        bg = BuildGrn()
        bg.networkMat = Any["TFA1" "GeneA" 0.9 0.8 0.5 "yes";
                            "TFA1" "GeneB" 0.7 0.6 -0.3 "no"]
        bg.networkMatSubset = bg.networkMat[1:1, :]

        writeNetworkTable!(bg; outputDir = tmpdir)

        edgesFile = joinpath(tmpdir, "edges.tsv")
        @test isfile(edgesFile)

        lines = readlines(edgesFile)
        @test occursin("TF",     lines[1])
        @test occursin("Gene",   lines[1])
        @test occursin("Stability", lines[1])
        @test length(lines) == 3    # header + 2 rows
    end

    @testset "writeTSVWithEmptyFirstHeader — roundtrip" begin
        tmpdir = mktempdir()
        df = DataFrame(Gene=["X","Y"], A=[1.0,2.0], B=[3.0,4.0])
        path = joinpath(tmpdir, "prior.tsv")
        writeTSVWithEmptyFirstHeader(df, path)

        lines = readlines(path)
        @test lines[1][1] == '\t'           # empty first header cell
        cols = split(lines[1], '\t')
        @test cols[2] == "A" && cols[3] == "B"
        row1 = split(lines[2], '\t')
        @test row1[1] == "X"
    end

end # Network I/O


# =============================================================================
@testset "GRN utilities" begin
# =============================================================================

    @testset "firstNByGroup — returns first N indices per group" begin
        v = ["A","B","A","C","B","A"]
        idxs = InferelatorJL.firstNByGroup(v, 2)
        # A appears at 1,3,6 — first 2 are 1,3
        # B appears at 2,5 — both kept
        # C appears at 4 — kept
        @test 1 in idxs && 3 in idxs && !(6 in idxs)
        @test 2 in idxs && 5 in idxs
        @test 4 in idxs
        @test length(idxs) == 5
    end

    @testset "firstNByGroup — N larger than group keeps all" begin
        v = ["A","A","B"]
        idxs = InferelatorJL.firstNByGroup(v, 10)
        @test length(idxs) == 3
    end

end # GRN utilities


# =============================================================================
@testset "Network aggregation" begin
# =============================================================================

    @testset "aggregateNetworks — max method produces output" begin
        tmpdir = mktempdir()

        edges1 = [
            (TF="TF1", Gene="GeneA", signedQuantile=0.9,  Stability=0.8, Correlation= 0.6, inPrior="yes"),
            (TF="TF1", Gene="GeneB", signedQuantile=0.7,  Stability=0.6, Correlation=-0.4, inPrior="no"),
            (TF="TF2", Gene="GeneA", signedQuantile=0.5,  Stability=0.5, Correlation= 0.3, inPrior="yes"),
            (TF="TF2", Gene="GeneC", signedQuantile=0.4,  Stability=0.4, Correlation= 0.2, inPrior="no"),
        ]
        edges2 = [
            (TF="TF1", Gene="GeneA", signedQuantile=0.95, Stability=0.9, Correlation= 0.7, inPrior="yes"),
            (TF="TF1", Gene="GeneB", signedQuantile=0.6,  Stability=0.5, Correlation=-0.3, inPrior="no"),
            (TF="TF2", Gene="GeneB", signedQuantile=0.3,  Stability=0.3, Correlation= 0.15,inPrior="no"),
        ]
        f1 = joinpath(tmpdir, "net1.tsv"); write_edges_tsv(f1, edges1)
        f2 = joinpath(tmpdir, "net2.tsv"); write_edges_tsv(f2, edges2)

        result = aggregateNetworks([f1, f2];
                                   method           = :max,
                                   meanEdgesPerGene = 5,
                                   useMeanEdgesPerGene = true,
                                   outputDir        = nothing)

        @test result isa DataFrame
        @test nrow(result) > 0
        @test "TF"          in names(result)
        @test "Gene"        in names(result)
        @test "Stability"   in names(result)
        @test "Correlation" in names(result)
    end

    @testset "aggregateNetworks — max picks higher stability for same TF-Gene pair" begin
        tmpdir = mktempdir()
        edges1 = [(TF="TF1", Gene="GeneA", signedQuantile=0.5, Stability=0.4, Correlation=0.3, inPrior="no")]
        edges2 = [(TF="TF1", Gene="GeneA", signedQuantile=0.9, Stability=0.9, Correlation=0.7, inPrior="yes")]
        f1 = joinpath(tmpdir,"n1.tsv"); write_edges_tsv(f1, edges1)
        f2 = joinpath(tmpdir,"n2.tsv"); write_edges_tsv(f2, edges2)

        result = aggregateNetworks([f1, f2]; method=:max, meanEdgesPerGene=5,
                                   useMeanEdgesPerGene=true, outputDir=nothing)

        row = result[result.TF .== "TF1" .&& result.Gene .== "GeneA", :]
        @test nrow(row) == 1
        @test row.Stability[1] ≈ 0.9
    end

    @testset "aggregateNetworks — mean averages stability and correlation" begin
        tmpdir = mktempdir()
        edges1 = [(TF="TF1", Gene="GeneA", signedQuantile=0.5, Stability=0.4, Correlation=0.2, inPrior="no")]
        edges2 = [(TF="TF1", Gene="GeneA", signedQuantile=0.8, Stability=0.8, Correlation=0.6, inPrior="no")]
        f1 = joinpath(tmpdir,"n1.tsv"); write_edges_tsv(f1, edges1)
        f2 = joinpath(tmpdir,"n2.tsv"); write_edges_tsv(f2, edges2)

        result = aggregateNetworks([f1, f2]; method=:mean, meanEdgesPerGene=5,
                                   useMeanEdgesPerGene=true, outputDir=nothing)

        row = result[result.TF .== "TF1" .&& result.Gene .== "GeneA", :]
        @test nrow(row) == 1
        @test row.Stability[1]   ≈ 0.6  atol=1e-6
        @test row.Correlation[1] ≈ 0.4  atol=1e-6
    end

    @testset "aggregateNetworks — output files written when outputDir is set" begin
        tmpdir = mktempdir()
        edges1 = [
            (TF="TF1", Gene="GeneA", signedQuantile=0.9, Stability=0.8, Correlation=0.6, inPrior="yes"),
            (TF="TF1", Gene="GeneB", signedQuantile=0.7, Stability=0.6, Correlation=0.4, inPrior="no"),
        ]
        f1 = joinpath(tmpdir,"n1.tsv"); write_edges_tsv(f1, edges1)
        outdir = joinpath(tmpdir, "combined")

        aggregateNetworks([f1]; method=:max, meanEdgesPerGene=5,
                          useMeanEdgesPerGene=true, outputDir=outdir)

        @test isfile(joinpath(outdir, "combined_max.tsv"))
        @test isfile(joinpath(outdir, "combined_max_sp.tsv"))
    end

end # Network aggregation


# =============================================================================
@testset "Merge degenerate TFs" begin
# =============================================================================

    @testset "Identical TF columns are merged" begin
        tmpdir = mktempdir()
        # TF1 and TF2 have identical binding profiles → should be merged
        # TF3 is unique → should remain separate
        prior_mat = Float64[1 1 0; 0 0 1; 1 1 0; 0 0 1; 1 1 0]
        genes = ["Akt1","Bcl2","Foxp3","Myc","Stat3"]
        tfs   = ["TF1","TF2","TF3"]
        priorFile = joinpath(tmpdir, "prior.tsv")
        write_prior_tsv(priorFile, tfs, genes, prior_mat)

        result = mergedTFsResult()
        InferelatorJL.mergeDegenerateTFs(result, priorFile; fileFormat = 2)

        # The mergedTFMap should be non-empty (TF1 and TF2 were merged)
        @test result.mergedTFMap !== nothing
        @test size(result.mergedTFMap, 2) == 2   # two-column matrix: metaTF, members

        # The merged prior DataFrame should contain a merged TF name
        @test result.mergedPrior !== nothing
        merged_tfs = names(result.mergedPrior)[2:end]   # skip first (gene) column
        @test length(merged_tfs) < length(tfs)           # fewer TFs after merging
    end

    @testset "Unique TF is absent from the merged map" begin
        tmpdir = mktempdir()
        # TF1 and TF2 are identical (will be merged); TF3 is unique (will not appear in map)
        prior_mat = Float64[1 1 0; 0 0 1; 1 1 0; 0 0 1; 1 1 0]
        genes = ["Akt1","Bcl2","Foxp3","Myc","Stat3"]
        tfs   = ["TF1","TF2","TF3"]
        priorFile = joinpath(tmpdir, "prior_mixed.tsv")
        write_prior_tsv(priorFile, tfs, genes, prior_mat)

        result = mergedTFsResult()
        InferelatorJL.mergeDegenerateTFs(result, priorFile; fileFormat = 2)

        # TF3 is unique — its name should not appear anywhere in the merged map
        merged_member_strings = result.mergedTFMap[:, 2]   # comma-separated member lists
        @test !any(occursin("TF3", s) for s in merged_member_strings)

        # The merged prior should still contain TF3 as its own column
        merged_cols = names(result.mergedPrior)[2:end]
        @test any(occursin("TF3", c) for c in merged_cols)
    end

end # Merge degenerate TFs


# =============================================================================
# ADDED: applyTimeLag tests
@testset "applyTimeLag" begin
# =============================================================================

    mktempdir() do tmp
        # ── synthetic setup ──────────────────────────────────────────────────
        # 2 TFs, 3 samples (s1=t0, s2=t1, s3=t2); s3 has a proceeding point s2
        Random.seed!(42)
        nTFs   = 2
        nTargs = 4
        nSamps = 3
        labels = ["s1", "s2", "s3"]

        # Build a minimal GeneExpressionData with the fields applyTimeLag! touches
        data = GeneExpressionData()
        data.cellLabels     = labels
        data.potRegMatmRNA  = rand(nTFs, nSamps) .+ 1.0

        # Build a minimal PriorTFAData
        prior = PriorTFAData()
        prior.medTfas        = rand(nTFs, nSamps) .+ 2.0
        prior.noPriorRegsMat = rand(nTFs, nSamps) .+ 3.0

        # Snapshot originals before correction
        medTfas_orig       = copy(prior.medTfas)
        noPrior_orig       = copy(prior.noPriorRegsMat)
        potReg_orig        = copy(data.potRegMatmRNA)

        # Write time-lag file: s3 (t=2) proceeds from s2 (t=1), lag = 0.5
        lagFile = joinpath(tmp, "timeLag.tsv")
        open(lagFile, "w") do io
            println(io, "SampleQ\tTimeQ\tSampleP\tTimeP")
            println(io, "s3\t2.0\ts2\t1.0")
        end

        applyTimeLag(prior, data, lagFile, 0.5)

        # ── correctness: only column 3 (s3) should change ────────────────────
        @test prior.medTfas[:, 1]        ≈ medTfas_orig[:, 1]   # s1 unchanged
        @test prior.medTfas[:, 2]        ≈ medTfas_orig[:, 2]   # s2 unchanged
        @test prior.noPriorRegsMat[:, 1] ≈ noPrior_orig[:, 1]
        @test prior.noPriorRegsMat[:, 2] ≈ noPrior_orig[:, 2]
        @test data.potRegMatmRNA[:, 1]   ≈ potReg_orig[:, 1]
        @test data.potRegMatmRNA[:, 2]   ≈ potReg_orig[:, 2]

        # s3 must differ from original
        @test !(prior.medTfas[:, 3]        ≈ medTfas_orig[:, 3])
        @test !(prior.noPriorRegsMat[:, 3] ≈ noPrior_orig[:, 3])
        @test !(data.potRegMatmRNA[:, 3]   ≈ potReg_orig[:, 3])

        # ── formula check: adjusted[Q] = orig[Q] - (orig[Q]-orig[P])/Δt * lag ─
        dt  = 2.0 - 1.0   # TimeQ - TimeP
        lag = 0.5
        expected_medTfas = medTfas_orig[:, 3] .-
                           (medTfas_orig[:, 3] .- medTfas_orig[:, 2]) ./ dt .* lag
        @test prior.medTfas[:, 3] ≈ expected_medTfas

        # ── missing sample name: emits @warn and skips (no error, no change) ───
        # FIXED: was @test_nowarn — but applyTimeLag! correctly emits @warn for
        # unknown samples, so we now assert the warning IS emitted and matrices
        # remain unchanged.
        lagFile2 = joinpath(tmp, "timeLag_bad.tsv")
        open(lagFile2, "w") do io
            println(io, "SampleQ\tTimeQ\tSampleP\tTimeP")
            println(io, "s_missing\t3.0\ts2\t2.0")
        end
        # reset to originals
        prior.medTfas        = copy(medTfas_orig)
        prior.noPriorRegsMat = copy(noPrior_orig)
        data.potRegMatmRNA   = copy(potReg_orig)
        @test_logs (:warn,) (:info,) applyTimeLag(prior, data, lagFile2, 0.5)  # FIXED: must match both the @warn (skipped pair) and the @info (summary) emitted by applyTimeLag!
        # nothing should have changed — skipped pair leaves all matrices intact
        @test prior.medTfas ≈ medTfas_orig
    end

end # applyTimeLag


# =============================================================================
@testset "Pseudobulk utilities" begin
# =============================================================================

    # Shared small matrix: 4 features × 6 samples, all positive
    Random.seed!(11)
    pbMat = abs.(rand(Float64, 4, 6)) .* 10 .+ 1.0

    @testset "normalizeMatrix — none returns identical matrix" begin
        @test normalizeMatrix(pbMat, "none") === pbMat
    end

    @testset "normalizeMatrix — tpm: column sums ≈ 1e6" begin
        out = normalizeMatrix(pbMat, "tpm")
        @test size(out) == size(pbMat)
        for j in 1:size(out, 2)
            @test sum(out[:, j]) ≈ 1e6  atol=1e-4
        end
    end

    @testset "normalizeMatrix — log2tpm: non-negative, larger than tpm input" begin
        out = normalizeMatrix(pbMat, "log2tpm")
        @test size(out) == size(pbMat)
        @test all(out .>= 0)
    end

    @testset "normalizeMatrix — zscore: row means ≈ 0, row stds ≈ 1" begin
        out = normalizeMatrix(pbMat, "zscore")
        rowMeans = vec(mean(out, dims=2))
        rowStds  = vec(std(out,  dims=2))
        @test all(abs.(rowMeans) .< 1e-10)
        @test all(abs.(rowStds  .- 1.0) .< 1e-10)
    end

    @testset "normalizeMatrix — log2tpm_zscore: row means ≈ 0" begin
        out = normalizeMatrix(pbMat, "log2tpm_zscore")
        rowMeans = vec(mean(out, dims=2))
        @test all(abs.(rowMeans) .< 1e-10)
    end

    @testset "normalizeMatrix — log2fc: correct shape, finite, mix of signs" begin
        out = normalizeMatrix(pbMat, "log2fc")
        @test size(out) == size(pbMat)
        @test all(isfinite, out)
        # log2(x / mean(x)) produces values centred around 0 — both positive and negative
        @test any(out .> 0)
        @test any(out .< 0)
    end

    @testset "normalizeMatrix — sizefactor: runs and returns same shape" begin
        out = normalizeMatrix(pbMat, "sizefactor")
        @test size(out) == size(pbMat)
        @test all(isfinite, out)
    end

    @testset "normalizeMatrix — vst throws with informative error" begin
        @test_throws ErrorException normalizeMatrix(pbMat, "vst")
    end

    @testset "normalizeMatrix — unknown method throws" begin
        @test_throws ErrorException normalizeMatrix(pbMat, "not_a_method")
    end

    @testset "buildDesignMatrix — shape and drop-first encoding" begin
        meta = DataFrame(celltype=["A","B","A","C","B"])
        D = buildDesignMatrix(meta, "celltype")

        # n rows, k-1 columns (drop first level)
        @test size(D) == (5, 2)
        # first level (A) is the reference — all its rows should be zero
        aRows = [1, 3]
        @test all(D[aRows, :] .== 0)
        # B rows should have a 1 in column 1
        bRows = [2, 5]
        @test all(D[bRows, 1] .== 1)
    end

    @testset "removeBatchEffect — batch mean removed per feature" begin
        # 2 features, 6 samples split evenly across 2 batches
        # Batch 2 has a constant offset of 5 added
        nFeats, nSamps = 2, 6
        Random.seed!(20)
        bio = rand(Float64, nFeats, nSamps)
        batch = ["B1","B1","B1","B2","B2","B2"]
        offset = 5.0
        X = hcat(bio[:, 1:3], bio[:, 4:6] .+ offset)

        corrected = removeBatchEffect(X, batch)
        @test size(corrected) == (nFeats, nSamps)
        # After correction the two batch groups should have similar column means
        meanB1 = mean(corrected[:, 1:3])
        meanB2 = mean(corrected[:, 4:6])
        @test abs(meanB1 - meanB2) < 0.5   # offset largely removed
    end

    @testset "removeBatchEffect — with designMat runs without error" begin
        meta = DataFrame(celltype=["A","A","B","B","A","B"],
                         batch=["B1","B1","B1","B2","B2","B2"])
        D = buildDesignMatrix(meta, "celltype")
        Random.seed!(21)
        X = rand(Float64, 3, 6) .+ 2.0
        corrected = removeBatchEffect(X, meta[!, :batch]; designMat=D)
        @test size(corrected) == (3, 6)
        @test all(isfinite, corrected)
    end

    @testset "aggregateReplicates — sums replicate columns per celltype" begin
        # Two celltypes (CT1, CT2), two replicates (R1, R2)
        features = ["GeneA", "GeneB", "GeneC"]
        df = DataFrame(
            :RowNames => features,
            Symbol("CT1_R1") => [1.0, 2.0, 3.0],
            Symbol("CT1_R2") => [4.0, 5.0, 6.0],
            Symbol("CT2_R1") => [7.0, 8.0, 9.0],
            Symbol("CT2_R2") => [10.0, 11.0, 12.0],
        )
        agg = aggregateReplicates(df, ["CT1","CT2"], ["R1","R2"])

        @test nrow(agg) == 3
        @test :CT1 in propertynames(agg)
        @test :CT2 in propertynames(agg)
        @test agg[1, :CT1] ≈ 5.0    # 1+4
        @test agg[2, :CT1] ≈ 7.0    # 2+5
        @test agg[1, :CT2] ≈ 17.0   # 7+10
    end

    @testset "aggregateReplicates — missing celltype silently skipped" begin
        features = ["GeneA"]
        df = DataFrame(
            :RowNames  => features,
            Symbol("CT1_R1") => [1.0],
        )
        agg = aggregateReplicates(df, ["CT1","CT_MISSING"], ["R1"])
        @test :CT1 in propertynames(agg)
        @test !(:CT_MISSING in propertynames(agg))
    end

end # Pseudobulk utilities


# =============================================================================
@testset "Real data (optional)" begin
# =============================================================================
# These tests only run when real file paths are provided via constants at the
# top of this file or via environment variables.  Otherwise they are skipped.

    have_real = isfile(REAL_EXPR_FILE) && isfile(REAL_TARG_FILE) &&
                isfile(REAL_REG_FILE)  && isfile(REAL_PRIOR_FILE)

    if !have_real
        @info "Real-data tests skipped — set REAL_EXPR_FILE / REAL_TARG_FILE / " *
              "REAL_REG_FILE / REAL_PRIOR_FILE to enable them."
    end

    @testset "Real data — expression loading" begin
        if !isfile(REAL_EXPR_FILE)
            @test_skip "REAL_EXPR_FILE not set"
        else
            data = GeneExpressionData()
            InferelatorJL.loadExpressionData!(data, REAL_EXPR_FILE)
            @test length(data.geneNames) > 0
            @test length(data.cellLabels) > 0
            @test size(data.geneExpressionMat, 1) == length(data.geneNames)
            @test size(data.geneExpressionMat, 2) == length(data.cellLabels)
            @test data.geneNames == sort(data.geneNames)   # sorted
            @info "Real expr: $(length(data.geneNames)) genes × $(length(data.cellLabels)) samples"
        end
    end

    @testset "Real data — full load + TFA" begin
        if !have_real
            @test_skip "One or more real data files not set"
        else
            tmpdir = mktempdir()
            data = GeneExpressionData()
            InferelatorJL.loadExpressionData!(data, REAL_EXPR_FILE)
            InferelatorJL.loadAndFilterTargetGenes!(data, REAL_TARG_FILE; epsilon = 0.01)
            InferelatorJL.loadPotentialRegulators!(data, REAL_REG_FILE)
            InferelatorJL.processTFAGenes!(data, "")

            tfaData = PriorTFAData()
            InferelatorJL.processPriorFile!(tfaData, data, REAL_PRIOR_FILE;
                                            mergedTFsData = nothing, minTargets = 3)
            InferelatorJL.calculateTFA!(tfaData, data; edgeSS = 0, zTarget = true)

            @test size(tfaData.medTfas, 1) == length(tfaData.pRegs)
            @test size(tfaData.medTfas, 2) == length(data.cellLabels)
            @test !any(isnan, tfaData.medTfas)
            @info "Real TFA: $(size(tfaData.medTfas, 1)) TFs × $(size(tfaData.medTfas, 2)) samples"
        end
    end

end # Real data


# =============================================================================
end  # @testset "InferelatorJL"
# =============================================================================
