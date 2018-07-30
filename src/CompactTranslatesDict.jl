__precompile__()
module CompactTranslatesDict

using BasisFunctions, CardinalBSplines, Domains

if VERSION < v"0.7-"
    nothing
else
    using LinearAlgebra
end

import Base: ==, getindex

using BasisFunctions: ShiftedIndex, ShiftedIndexList
import BasisFunctions: length, is_biorthogonal, is_basis, name, ordering
import BasisFunctions: has_unitary_transform, support, has_grid, grid, period
import BasisFunctions: stepsize, has_grid_transform, compatible_grid, approx_length
import BasisFunctions: native_nodes, transform_from_grid, transform_to_grid
import BasisFunctions: grid_evaluation_operator, unsafe_eval_element
import BasisFunctions: Gram, UnNormalizedGram, dual, dict_in_support
import BasisFunctions: matrix, apply!, adjoint

export is_biorthogonal, is_basis, ordering, has_unitary_transform, support, has_grid
export grid, period, Gram

export CompactPeriodicTranslationDict, dual, discrete_dual
export BSplineTranslatesBasis, bspline_platform
export IndexableHorizontalBandedOperator, IndexableVerticalBandedOperator, matrix

include("operators/banded_operators.jl")
include("operators/ext_res_operator.jl")
include("translates.jl")
include("translates_of_bsplines.jl")
include("approximation.jl")
end
