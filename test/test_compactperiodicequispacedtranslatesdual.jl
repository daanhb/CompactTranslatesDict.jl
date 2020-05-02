

using Test, BasisFunctions, DomainSets, CompactTranslatesDict, CardinalBSplines
using GridArrays: similargrid

N = 10
x = LinRange(-3,3,100)
B = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(N,-1.123,1.432), x->CenteredBSpline(2)(N*x),(-2..2)/N)

m = 2
D = CompactPeriodicEquispacedTranslatesDual(B, m)
m1 = evaluation(B,similargrid(interpolation_grid(B), Float64,m*N))
m2 = evaluation(D,similargrid(interpolation_grid(B), Float64,m*N))
@test m1'm2≈IdentityOperator(B, B)
μ = discretemeasure(similargrid(interpolation_grid(B),Float64,m*N))
@test mixedgram(B,D,μ)≈IdentityOperator(B, B)

using InfiniteVectors
b = CompactPeriodicEquispacedTranslatesDuals.signal(B,2)
primal_signal = PeriodicInfiniteVector(b, 20)[0:19]
c = inv(b, 2,K=D.minimalK)'
dual_signal = PeriodicInfiniteVector(c, 20)[0:19]

@test evaluation(B, points(μ)).A[:,1]≈primal_signal
@test evaluation(D, points(μ)).A[:,1]≈dual_signal
