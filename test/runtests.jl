# test/runtests.jl
using Test
using InferelatorJL

@testset "InferelatorJL" begin

    @testset "Structs instantiate" begin
        @test GeneExpressionData() isa GeneExpressionData
        @test PriorTFAData()       isa PriorTFAData
        @test mergedTFsResult()    isa mergedTFsResult
        @test GrnData()            isa GrnData
        @test BuildGrn()           isa BuildGrn
    end

    # TODO: add data-loading tests once example data is available
    # TODO: add TFA estimation tests
    # TODO: add GRN regression tests

end
