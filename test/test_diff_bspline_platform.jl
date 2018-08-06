

using CompactTranslatesDict, Base.Test, CompactTranslatesDict.SymbolicDifferentialOperators

CompactTranslatesDict.bspline_basis(Float64, (10,10), (2,3), (δx^2+δy^3)*(δx^2+δy^3))


B = BSplineTranslatesBasis(10,2)⊗BSplineTranslatesBasis(10,3)
L = CompactTranslatesDict.bspline_evaluation_operator(Float64, (10,10), (2,3), BasisFunctions.grid(B),(δx^2+δy^3)*(δx^2+δy^3))

@test CompactTranslatesDict.is_fastly_indexable(L)
