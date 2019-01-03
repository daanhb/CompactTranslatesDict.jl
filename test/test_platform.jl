
using CompactTranslatesDict, CompactTranslatesDict.SymbolicDifferentialOperators, BasisFunctions, StaticArrays, DomainSets
types = [Float64, BigFloat]
(VERSION<v"0.7-") ? using Base.Test : using Test

function test_bspline_platform(T)
    #  1D
    init = 4
    for oversampling in [1,2,4],  degree in 2:3, i in [1,3]
        platform = bspline_platform(T, init, degree, oversampling)

        P = primal(platform,i)
        D = dual(platform,i)
        S = BasisFunctions.sampler(platform,i)

        B = P
        g = BasisFunctions.oversampled_grid(B, oversampling)
        E = CirculantOperator(evaluation_matrix(B[1],g)[:])*IndexExtensionOperator(B,gridbasis(g),1:oversampling:length(g))
        G = CirculantOperator(E'E*[1,zeros(length(g)-1)...]/length(B))
        DG = BasisFunctions.wrap_operator(B, B, inv(G))

        e = map(T,rand(length(B)))
        @test evaluation_operator(D,g)*e≈evaluation_matrix(D,g)*e
        @test evaluation_operator(D,g)*e≈evaluation_operator(P,g)*(matrix(DG)*e)
        @test evaluation_operator(D,g)'*evaluation_operator(P,g)*e ≈length(B)e
        @test S*exp≈broadcast(exp,g)
    end

    #  ND
    init = (3,4)
    degree = (1,3)
    T = Float64
    oversampling = 2
    for oversampling in [1,2,4], i in [1,2]
        platform = bspline_platform(T, init, degree, oversampling)

        P = primal(platform,i)
        D = dual(platform,i)
        S = BasisFunctions.sampler(platform,i)


        B = P
        B1, B2 = elements(P)
        g1 = BasisFunctions.oversampled_grid(B1,oversampling)
        g2 = BasisFunctions.oversampled_grid(B2,oversampling)
        g = g1×g2

        E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(B1,gridbasis(g1),1:oversampling:length(g1))
        E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(B2,gridbasis(g2),1:oversampling:length(g2))

        G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]/length(B1));
        G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]/length(B2));
        G = G1⊗G2

        DG = BasisFunctions.wrap_operator(B, B, inv(G))

        e = map(T,rand(size(B)...))
        @test (evaluation_operator(D,g)*e)[:]≈evaluation_matrix(D,g)*e[:]
        @test (evaluation_operator(D,g)*e)[:]≈(evaluation_operator(P,g)*reshape(matrix(DG)*e[:],size(P)...))[:]
        @test evaluation_operator(D,g)'*evaluation_operator(P,g)*e ≈length(B)e
        f = (x,y)->exp(x*y)
        @test S*f≈broadcast(f,g)
    end

    init = (3,3)
    degree = (1,3)
    oversampling = 2
    center = @SVector [.5,.5]
    domain  = disk(.3,center)
    platform=bspline_platform(T, init, degree, oversampling)
    i = 2
    B = primal(platform,i)
    D = dual(platform,i)
    g = grid(sampler(platform,i))
    Aop = BasisFunctions.A(platform, i)
    Zop = BasisFunctions.Zt(platform, i)'


    e = map(T,rand(size(B)...))
    @test evaluation_operator(B, g)*e≈Aop*e
    @test evaluation_operator(D, g)*e≈Zop*e*length(D)
    @test BasisFunctions.Zt(platform, i)*Aop*e ≈ e
end

function test_diff_bspline_platform(T)

    pg = CompactTranslatesDict.primal_diff_bspline_generator(Float64, (3,3), δx^2*δy + 1, ('x','y'))
    P = CompactTranslatesDict.bspline_basis(Float64, (9,10), (3,3), δx^2*δy + 1, ('x','y'))
    dg = CompactTranslatesDict.dual_diff_bspline_generator(pg, 2)
    D = dg((9,10))

    e = rand(P)
    @test evaluation_operator(pg((9,10)))*e≈evaluation_operator(P)*e

    ED = evaluation_operator(D, oversampling=2)
    EP = evaluation_operator(P, oversampling=2)
    G = DiscreteGram(P, oversampling=2)
    @test EP*inv(G)*e≈ED*e

    platform = bspline_platform(T, (4,4), (3,3), 2, δx*δy + 1)
    primal(platform, 1)
    dual(platform, 1)
    sampler(platform, 1)
    a = A(platform, 1)
    zt = Zt(platform, 1)

    e = rand(src(a))
    @test zt*a*e≈e
end


for T in types
    @testset "$(rpad("Translates of B spline platforms",80))" begin
        test_bspline_platform(T)
    end
end

T = Float64
@testset "$(rpad("Translates of B spline platforms on differential equation",80))" begin
    test_diff_bspline_platform(T)
end
