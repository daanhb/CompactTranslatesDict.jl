using CompactTranslatesDict

include("test_compact_approximation.jl")
include("test_bsplinetranslatedbasis.jl")

for T in [Float64, BigFloat]
    @testset "$(rpad("Periodic translate expansions",80))" begin
        test_generic_periodicbsplinebasis(T) end

    @testset "$(rpad("Translates of B spline expansions",80))" begin
        test_translatedbsplines(T)
        test_translatedsymmetricbsplines(T)
        # test_orthonormalsplinebasis(T)
        # test_discrete_orthonormalsplinebasis(T)
        test_dualsplinebasis(T)
        test_discrete_dualsplinebasis(T)
        test_bspline_platform(T)
        test_sparsity_speed(T)
    end
end
