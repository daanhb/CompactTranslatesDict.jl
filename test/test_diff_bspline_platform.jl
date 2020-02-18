

using CompactTranslatesDict, BasisFunctions, CompactTranslatesDict.SymbolicDifferentialOperators
using Test
CompactTranslatesDict.bspline_basis(Float64, (10,10), (2,3), (δx^2+δy^3)*(δx^2+δy^3))


B = BSplineTranslatesBasis(10,2)⊗BSplineTranslatesBasis(10,3)
L = evaluation(CompactTranslatesDict.bspline_basis(Float64, (10,10), (2,3), (δx^2+δy^3)*(δx^2+δy^3)))

@test CompactTranslatesDict.is_fastly_indexable(L)
