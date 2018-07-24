__precompile__()
module CompactTranslatesDict

using BasisFunctions, CardinalBSplines, Domains

if VERSION < v"0.7-"
    nothing
else
    using LinearAlgebra
end

import Base.==

using BasisFunctions: ShiftedIndex, ShiftedIndexList
import BasisFunctions: length, is_biorthogonal, is_basis, name, ordering
import BasisFunctions: has_unitary_transform, support, has_grid, grid, period
import BasisFunctions: stepsize, has_grid_transform, compatible_grid, approx_length
import BasisFunctions: native_nodes, transform_from_grid, transform_to_grid
import BasisFunctions: grid_evaluation_operator, unsafe_eval_element
import BasisFunctions: Gram, UnNormalizedGram, dual, dict_in_support

export is_biorthogonal, is_basis, ordering, has_unitary_transform, support, has_grid
export grid, period, Gram

# from bases/translates/translation_dict.jl
export CompactPeriodicTranslationDict, dual, discrete_dual
# from bases/translates/translates_of_bsplines.jl
export BSplineTranslatesBasis, bspline_platform

include("translates.jl")
include("translates_of_bsplines.jl")
include("approximation.jl")
end
