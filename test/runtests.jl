using CompactTranslatesDict
types = [Float64, BigFloat]





using CompactTranslatesDict, Test, DomainSets
@testset begin
    g1 = GenericTranslates(EquispacedGrid(10,0,1), exp)
    g2 = GenericEquispacedTranslates(PeriodicEquispacedGrid(10,0,1), exp)
    g3 = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(10,0,1), exp)
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

@testset begin
    test_bspline_orthogonality_orthonormality()
end

for T in types
    @testset "$(rpad("Translates of B spline expansions",80))" begin
        test_generic_periodicbsplinebasis(BSplineTranslatesBasis, T)
        test_translatedbsplines(T)
    end
end

@testset "gramoperator" begin
    B = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(10,0,2π), cos)
    μ = grid(discretemeasure(PeriodicEquispacedGrid(10,0,2π)))
    g1 = gramoperator(B, discretemeasure(PeriodicEquispacedGrid(10,0,2π)))
    g2 = BasisFunctions.default_gramoperator(B, discretemeasure(PeriodicEquispacedGrid(10,0,2π)))
    @test g1≈g2

    B = BSplineTranslatesBasis(10,3)
    g0 = BasisFunctions.default_gramoperator(B;overquad=100)
    g1 = gramoperator(B;overquad=100)
    g2 = gramoperator(B, FourierMeasure();overquad=100)
    g3 = BasisFunctions.default_gramoperator(B, FourierMeasure();overquad=100)
    @test g1≈g0≈g2≈g3
    @test g1 isa CirculantOperator
    @test g2 isa CirculantOperator

    g0 = gramoperator(B, discretemeasure(BasisFunctions.interpolation_grid(B)))
    g1 = BasisFunctions.default_gramoperator(B, discretemeasure(BasisFunctions.interpolation_grid(B)))
    @test g0≈g1
    @test g0 isa CirculantOperator
    @test BasisFunctions.hasinterpolationgrid(B) ≈ hasinterpolationgrid(B)
end

# include("test_platform.jl")
