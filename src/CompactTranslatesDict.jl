__precompile__()
module CompactTranslatesDict


include("SymbolicDifferentialOperators/SymbolicDifferentialOperators.jl")
using .SymbolicDifferentialOperators: SumDifferentialOperator, ProductDifferentialOperator, PartialDifferentialOperator,IdentityDifferentialOperator, AbstractDiffOperator
using .SymbolicDifferentialOperators: dimension_names, coefficient, operator, scalar, dimension_name

using BasisFunctions, CardinalBSplines, Domains

if VERSION < v"0.7-"
    nothing
else
    using LinearAlgebra
end

import Base: ==, getindex

using BasisFunctions: ShiftedIndex, ShiftedIndexList, eigenvalues, product_domaintype, promote_coeftype, ProductGrid
import BasisFunctions: length, is_biorthogonal, is_basis, name, ordering
import BasisFunctions: has_unitary_transform, support, has_grid, grid, period
import BasisFunctions: stepsize, has_grid_transform, compatible_grid, approx_length
import BasisFunctions: native_nodes, transform_from_grid, transform_to_grid
import BasisFunctions: grid_evaluation_operator, unsafe_eval_element
import BasisFunctions: Gram, UnNormalizedGram, dual, dict_in_support
import BasisFunctions: matrix, apply!, adjoint, unsafe_wrap_operator
import BasisFunctions: instantiate, extension_operator, restriction_operator, resize
import BasisFunctions: elements, native_index, src, dest

export is_biorthogonal, is_basis, ordering, has_unitary_transform, support, has_grid
export grid, period, Gram

export CompactPeriodicTranslationDict, dual, discrete_dual
export BSplineTranslatesBasis, bspline_platform, DiffBSplineTranslatesBasis
export degree, Bdiff
export IndexableHorizontalBandedOperator, IndexableVerticalBandedOperator, matrix
export CompactTranslatesTensorProductDict, BSplineTensorProductDict, evaluation_operator

export ⊗


include("operators/banded_operators.jl")
include("operators/sumoperator.jl")
include("operators/ext_res_operator.jl")
include("dictionary/translates.jl")
include("dictionary/translates_of_bsplines.jl")

include("dictionary/tensor.jl")
include("dictionary/sum.jl")
include("platform/bspline_platform.jl")
include("platform/diff_bspline_platform.jl")
end
