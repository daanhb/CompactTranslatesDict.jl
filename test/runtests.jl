using CompactTranslatesDict

if (VERSION < v"0.7-")
    types = [Float64, BigFloat]
else
    types = (Float64,)
end
include("test_bsplinetranslatedbasis.jl")

for T in types
    @testset "$(rpad("Translates of B spline expansions",80))" begin
        test_generic_periodicbsplinebasis(BSplineTranslatesBasis, T)
        test_translatedbsplines(T)
        test_bspline_platform(T)
        test_sparsity_speed(T)
    end
end
