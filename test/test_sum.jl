
using CompactTranslatesDict, DomainSets, BasisFunctions, Test

@testset "$(rpad("Summation of tensorproduct",80))" begin
    B = BSplineTranslatesBasis(10,2)⊗DiffBSplineTranslatesBasis(10,3,1)

    @test isa(B, CompactTranslatesTensorProductDict)
    @test isa(B, BSplineTensorProductDict)

    @test BSplineTensorProductDict <: CompactTranslatesTensorProductDict

    @test degree(B) == (2,3)

    @test Bdiff(B) == (0,1)

    BB = CompactTranslatesDict.CompactTranslationDictSum((B, B, B), [.2,.4,.7])
    S = CompactTranslatesDict.grid_evaluation_operator(BB, GridBasis(BB), interpolation_grid(BB))
    e = rand(src(S))
    @test S*e ≈ 1.3CompactTranslatesDict.grid_evaluation_operator(B, GridBasis(B), interpolation_grid(B))*e

    @test  1.3evaluation_operator(B, interpolation_grid(B))*e ≈ evaluation_operator(BB, interpolation_grid(BB))*e
    @test matrix(S)≈matrix(S[:,:])

end
