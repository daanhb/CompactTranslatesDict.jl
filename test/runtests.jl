using CompactTranslatesDict
types = [Float64, BigFloat]

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

# include("test_platform.jl")
