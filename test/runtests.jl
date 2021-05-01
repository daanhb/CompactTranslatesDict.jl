types = (Float64,)

using CompactTranslatesDict, BasisFunctions, Test, DomainSets
using BasisFunctions: period, isperiodic

@testset "generic translates" begin
    g1 = GenericTranslates(EquispacedGrid(10,0,1), exp)
    g2 = GenericEquispacedTranslates(PeriodicEquispacedGrid(10,0,1), exp)
    g3 = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(10,0,1), exp, Interval(0,.5))
    for g in (g1,g2,g3)
        @test support(g) ≈ UnitInterval()
    end
    for g in (g2,g3)
        @test step(g) ≈ 1/10
    end
    @test isperiodic(g3)
    @test period(g3) ≈ 1
end



include("test_bsplinetranslatedbasis.jl")

@testset "orthogonality and orthonormality" begin
    test_bspline_orthogonality_orthonormality()
end

for T in types
    @testset "Translates of B spline expansions" begin
        test_generic_periodicbsplinebasis(BSplineTranslatesBasis, T)
        test_translatedbsplines(T)
    end
end

@testset "gram" begin
    B = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(10,0,2π), cos, Interval(0.,.4))
    μ = discretemeasure(PeriodicEquispacedGrid(10,0,2π))
    g1 = gram(B, μ)
    g2 = BasisFunctions.default_gram(B, μ)
    @test g1≈g2

    B = BSplineTranslatesBasis(10,3)
    g0 = BasisFunctions.default_gram(B;overquad=100)
    g1 = gram(B;overquad=100)
    g2 = gram(B, FourierWeight();overquad=100)
    g3 = BasisFunctions.default_gram(B, FourierWeight();overquad=100)
    @test g1≈g0≈g2≈g3
    @test g1 isa CirculantOperator
    @test g2 isa CirculantOperator

    g0 = gram(B, discretemeasure(BasisFunctions.interpolation_grid(B)))
    g1 = BasisFunctions.default_gram(B, discretemeasure(BasisFunctions.interpolation_grid(B)))
    @test g0≈g1
    @test g0 isa CirculantOperator
    @test BasisFunctions.hasinterpolationgrid(B) ≈ hasinterpolationgrid(B)
end

using Test
@testset "Compact Duals" begin
    include("test_compactperiodicequispacedtranslatesdual.jl")
end
