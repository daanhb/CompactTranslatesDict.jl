module CompactTranslatesDict

using BasisFunctions, CardinalBSplines, Domains
import Base.==

# from bases/translates/translation_dict.jl
export CompactPeriodicTranslationDict, dual, discrete_dual
# from bases/translates/translates_of_bsplines.jl
export BSplineTranslatesBasis, SymBSplineTranslatesBasis, OrthonormalSplineBasis, DiscreteOrthonormalSplineBasis
export bspline_platform

include("translates.jl")
include("translates_of_bsplines.jl")
include("compact_approximation.jl")
end
