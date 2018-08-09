

using CompactTranslatesDict, Base.Test, Domains, BasisFunctions

@testset "$(rpad("Summation of tensorproduct",80))" begin
    B = BSplineTranslatesBasis(10,2)⊗DiffBSplineTranslatesBasis(10,3,1)

    @test isa(B, CompactTranslatesTensorProductDict)
    @test isa(B, BSplineTensorProductDict)

    @test BSplineTensorProductDict <: CompactTranslatesTensorProductDict

    @test degree(B) == (2,3)

    @test Bdiff(B) == (0,1)

    BB = CompactTranslatesDict.CompactTranslationDictSum((B, B, B), [.2,.4,.7])
    S = CompactTranslatesDict.grid_evaluation_operator(BB, BasisFunctions.gridbasis(BB), BasisFunctions.grid(BB))
    e = rand(BasisFunctions.src(S))
    @test S*e ≈ 1.3CompactTranslatesDict.grid_evaluation_operator(B, BasisFunctions.gridbasis(B), BasisFunctions.grid(B))*e

    @test 1.3evaluation_operator(B)*e ≈ evaluation_operator(BB)*e
    @test matrix(S)≈matrix(S[:,:])

end
