

using CompactTranslatesDict, Base.Test, SymbolicDifferentialEquations

CompactTranslatesDict.generate_bspline_basis(Float64, (10,10), (2,3), (δx^2+δy^3)*(δx^2+δy^3))
