using BasisFunctions, BasisFunctions.Test
using CompactTranslatesDict

if VERSION < v"0.7-"
    using Base.Test
else
    using Test
end
@testset "$(rpad("test banded operators",80))" begin
    ELT = Float64
    a = [ELT(1),ELT(2),ELT(3)]
    H = IndexableHorizontalBandedOperator(FourierBasis{ELT}(6,0,1), FourierBasis{ELT}(3,0,1),a,3,2)
    test_generic_operator_interface(H, ELT)
    h = matrix(H)
    e = zeros(ELT,6,1)
    for i in 1:3
        e .= 0
        I = (2+(i-1)*3) .+ (1:3)
        for (j,k) in enumerate(mod.(I .- 1,6) .+ 1)
            e[k] = a[j]
        end
        @test e ≈ h[i,:]
    end
    matrix(H')≈matrix(H)'

    V = IndexableVerticalBandedOperator(FourierBasis{ELT}(3,0,1), FourierBasis{ELT}(6,0,1),a,3,2)
    test_generic_operator_interface(V, ELT)
    v = matrix(V)
    e = zeros(ELT,6,1)
    for i in 1:3
        e .= 0
        I = (2+(i-1)*3) .+ (1:3)
        for (j,k) in enumerate(mod.(I .- 1,6) .+ 1)
            e[k] = a[j]
        end
        @test e ≈ v[:,i]
    end

    matrix(V')≈matrix(V)'
end


@testset "$(rpad("test indexable operators",80))" begin
    d = BSplineTranslatesBasis(10,3)
    s = BSplineTranslatesBasis(20,3)

    M = IndexableVerticalBandedOperator(d, s, [1,2], 2, 1)
    C = CompactTranslatesDict.ExtResOperator(1:3, M, 1:5)
    c = matrix(C)

    @test c == matrix(M)[1:3, 1:5]
    @test matrix(M[1:3, 2:4]) == matrix(M)[1:3,2:4]
    @test matrix(M[:,1:3]) == matrix(M)[:,1:3]
    @test matrix(M[2:3,:]) == matrix(M)[2:3,:]
    @test c == [C[i,j] for i in 1:size(C, 1), j in 1:size(C,2)]
    @test matrix(C[1:2,1:2]) == c[1:2,1:2]
    @test matrix(C[:,1:2])==c[:,1:2]
    @test matrix(C[1:2,:])==c[1:2,:]


    M = evaluation_operator(s⊗d,oversampling=2)
    i = [CartesianIndex(1,1), CartesianIndex(2,3), CartesianIndex(3,1)]
    j = [CartesianIndex(1,1), CartesianIndex(2,3), CartesianIndex(3,1), CartesianIndex(5,6), CartesianIndex(1,10)]

    k = linear_index.(dest(M), i)
    l = linear_index.(src(M), j)

    C = CompactTranslatesDict.ExtResOperator(i, M, j)
    c = matrix(C)


    @test c == matrix(M)[k,l]
    @test matrix(M[i, j]) == matrix(M)[k,l]
    @test matrix(M[:,j]) == matrix(M)[:,l]
    @test matrix(M[i,:]) == matrix(M)[k,:]
    @test c == [C[i,j] for i in 1:size(C, 1), j in 1:size(C,2)]
    @test matrix(C[1:2,1:2]) == c[1:2,1:2]
    @test matrix(C[:,1:2])==c[:,1:2]
    @test matrix(C[1:2,:])==c[1:2,:]
end
