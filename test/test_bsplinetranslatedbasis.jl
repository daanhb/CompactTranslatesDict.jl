using BasisFunctions, DomainSets, CompactTranslatesDict, Test, LinearAlgebra
using BasisFunctions.Test: test_orthogonality_orthonormality

function test_generic_periodicbsplinebasis(B,T)
    tol = sqrt(eps(real(T)))
    n = 5
    b = B(n,3, T)
    @test leftendpoint(support(b)) == 0
    @test rightendpoint(support(b))==1

    @test length(b)==5
    @test CompactTranslatesDict.degree(b)==3
    @test is_basis(b)
    @test BasisFunctions.is_biorthogonal(b)
    @test !BasisFunctions.is_orthogonal(b)
    @test !BasisFunctions.is_orthonormal(b)
    @test !has_unitary_transform(b)

    @test infimum(CompactTranslatesDict.kernel_span(b)) <= 0 <= supremum(CompactTranslatesDict.kernel_span(b))
    @test 0 < supremum(CompactTranslatesDict.kernel_span(b)) - infimum(CompactTranslatesDict.kernel_span(b)) < BasisFunctions.period(b)

    @test instantiate(B, 4, Float16)==B(4,3,Float16)
    @test resize(b, 20)==B(20,CompactTranslatesDict.degree(b),T)
    @test interpolation_grid(b)==PeriodicEquispacedGrid(n,0,1)
    @test BasisFunctions.period(b)==T(1)
    @test BasisFunctions.stepsize(b)==T(1//5)
end

function test_translatedbsplines(T)

    tol = sqrt(eps(real(T)))
    n = 5
    bb = BSplineTranslatesBasis(n, 1, T; scaled=true)
    b = BSplineTranslatesBasis(n, 1, T; scaled=false)
    e = rand(n)
    @test norm(Gram(b)*e-Gram(bb)*e/n) < tol

    b = BSplineTranslatesBasis(n,3, T)

    @test BasisFunctions.name(b) == "Dictionary of equispaced translates of a kernel function (B spline of degree 3)"

    @test infimum(support(b,1))≈ 0
    @test infimum(support(b,2))≈1//5
    @test infimum(support(b,3))≈0
    @test infimum(support(b,4))≈0
    @test infimum(support(b,5))≈0
    @test supremum(support(b,1))≈4//5
    @test supremum(support(b,2))≈1
    @test supremum(support(b,3))≈1
    @test supremum(support(b,4))≈1
    @test supremum(support(b,5))≈1

    t = .001
    @test in_support(b,1,.5)
    @test !in_support(b,1,.8+t)
    @test !in_support(b,1,1. -t)
    @test in_support(b,3,.2-t)
    @test in_support(b,3,.4+t)
    @test !in_support(b,3,.2+t)
    @test !in_support(b,3,.4-t)

    @test BasisFunctions.compatible_grid(b, grid(b))
    @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1))
    @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n+1,0,1))
    @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0.1,1))
    @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1.1))

    grid(BSplineTranslatesBasis(n,2, T)) == MidpointEquispacedGrid(n,0,1)
    @test CompactTranslatesDict.degree(BSplineTranslatesBasis(5,2, T)) == 2
    b = BSplineTranslatesBasis(n,2,T)
    @test BasisFunctions.compatible_grid(b, grid(b))
    @test !BasisFunctions.compatible_grid(b, PeriodicEquispacedGrid(n,0,1))
    @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n+1,0,1))
    @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0.1,1))
    @test !BasisFunctions.compatible_grid(b, MidpointEquispacedGrid(n,0,1.1))

    # Test extension_operator and invertability of restriction_operator w.r.t. extension_operator.
    n = 8
    for degree in 0:3
        b = BSplineTranslatesBasis(n, degree, T)
        basis_ext = extend(b)
        r = restriction_operator(basis_ext, b)
        e = extension_operator(b, basis_ext)
        @test (VERSION < v"0.7-") ? abs(sum(eye(n)-matrix(r*e))) < tol : abs(sum(Matrix(1.0I, n, n) -matrix(r*e))) < tol

        grid_ext = grid(basis_ext)
        L = evaluation_operator(b, grid_ext)
        e = random_expansion(b)
        z = L*e
        L2 = evaluation_operator(basis_ext, grid_ext) * extension_operator(b, basis_ext)
        z2 = L2*e
        @test maximum(abs.(z-z2)) < tol
    end

    println("Expect 16 warnings")
    if T == Float64
        for K in 0:3
            for s2 in 5:6
                s1 = s2<<1
                b1 = BSplineTranslatesBasis(s1,K,T)
                b2 = BSplineTranslatesBasis(s2,K,T)

                e1 = random_expansion(b1)
                e2 = random_expansion(b2)

                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b2, GridBasis(b2))*e2) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b2, GridBasis(b2), grid(b2))*e2)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b2, GridBasis(b1))*e2) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b2, GridBasis(b1), grid(b1))*e2)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b1, GridBasis(b1))*e1) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b1, GridBasis(b1), grid(b1))*e1)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b1, GridBasis(b2))*e1) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b1, GridBasis(b2), grid(b2))*e1)

                mr = matrix(restriction_operator(b1, b2))
                me = matrix(extension_operator(b2, b1))
                pinvme = pinv(me)
                r = rand(size(pinvme,2))
                @test pinvme*r ≈ mr*r
            end
        end

        for K in 0:3
            for s2 in 5:6
                s1 = s2<<1
                b1 = BSplineTranslatesBasis(s1,K,T; scaled=true)
                b2 = BSplineTranslatesBasis(s2,K,T; scaled=true)

                e1 = random_expansion(b1)
                e2 = random_expansion(b2)

                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b2, GridBasis(b2))*e2) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b2, GridBasis(b2), grid(b2))*e2)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b2, GridBasis(b1))*e2) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b2, GridBasis(b1), grid(b1))*e2)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b1, GridBasis(b1))*e1) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b1, GridBasis(b1), grid(b1))*e1)
                @test BasisFunctions.coefficients(BasisFunctions.default_evaluation_operator(b1, GridBasis(b2))*e1) ≈ BasisFunctions.coefficients(grid_evaluation_operator(b1, GridBasis(b2), grid(b2))*e1)

                mr = matrix(restriction_operator(b1, b2))
                me = matrix(extension_operator(b2, b1))
                pinvme = pinv(me)
                r = rand(size(pinvme,2))
                @test pinvme*r ≈ mr*r
            end
        end
    end

    @test_throws AssertionError restriction_operator(BSplineTranslatesBasis(4,0,T), BSplineTranslatesBasis(3,0,T))
    @test_throws AssertionError extension_operator(BSplineTranslatesBasis(4,0,T), BSplineTranslatesBasis(6,0,T))
end


function test_bspline_orthogonality_orthonormality()
    for m in [FourierMeasure(),
                BasisFunctions.DiscreteMeasure(PeriodicEquispacedGrid(4,0,1)),
                BasisFunctions.DiscreteMeasure(MidpointEquispacedGrid(4,0,1)),
                BasisFunctions.DiscreteMeasure(PeriodicEquispacedGrid(8,0,1)),
                BasisFunctions.DiscreteMeasure(MidpointEquispacedGrid(8,0,1))]
        @test BasisFunctions.unsafe_matrix(gramoperator(B, m)) isa Circulant
        test_orthogonality_orthonormality(B, false, false, m)
    end
end



# using Plots, BasisFunctions, CompactTranslatesDict
# n = 7
# k = 0
# b = BSplineTranslatesBasis(n,k)
# t = linspace(-0,2,200)
# f = map(x->BasisFunctions.unsafe_eval_element(b,1,x), t)
#
# plot!(t,f)
# f = CompactTranslatesDict.eval_kernel.(b, t)
# plot!(t,f)
# f = map(x->BasisFunctions.eval_expansion(b,ones(n),x),t)
# # f = map(x->BasisFunctions.eval_expansion(b,[1,0,0,0,0,0,0,0,0,0],x),t)
# plot!(t,f,ylims=[-n,2n])
