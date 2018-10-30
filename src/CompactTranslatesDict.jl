__precompile__()
module CompactTranslatesDict

include("SymbolicDifferentialOperators/SymbolicDifferentialOperators.jl")
using .SymbolicDifferentialOperators: SumDifferentialOperator, ProductDifferentialOperator, PartialDifferentialOperator,IdentityDifferentialOperator, AbstractDiffOperator, ScaledDifferentialOperator
using .SymbolicDifferentialOperators: dimension_names, coefficient, operator, scalar, dimension_name

using BasisFunctions, CardinalBSplines, Domains

if VERSION < v"0.7-"
    using Compat
    using Compat: @warn
else
    using LinearAlgebra, FFTW
end

import Base: ==, getindex

using BasisFunctions: ShiftedIndex, ShiftedIndexList, eigenvalues, product_domaintype, promote_coeftype, ProductGrid
using BasisFunctions: forward_fourier_operator, op_eltypes, DictionaryOperator, ×
import BasisFunctions: length, is_biorthogonal, is_basis, name, ordering, oversampled_grid
import BasisFunctions: has_unitary_transform, support, has_grid, grid, period
import BasisFunctions: stepsize, has_grid_transform, compatible_grid, approx_length
import BasisFunctions: native_nodes, transform_from_grid, transform_to_grid
import BasisFunctions: grid_evaluation_operator, unsafe_eval_element, similar_dictionary
import BasisFunctions: Gram, UnNormalizedGram, dual, dict_in_support
import BasisFunctions: matrix, apply!, adjoint, unsafe_wrap_operator, matrix!
import BasisFunctions: instantiate, extension_operator, restriction_operator, resize
import BasisFunctions: elements, native_index, src, dest, *, +, -, inv, pinv, element

export is_biorthogonal, is_basis, ordering, has_unitary_transform, support, has_grid
export grid, period, Gram

export CompactPeriodicTranslationDict, dual, discrete_dual
export BSplineTranslatesBasis, DiffBSplineTranslatesBasis
export degree, Bdiff, src, dest
export IndexableHorizontalBandedOperator, IndexableVerticalBandedOperator, matrix
export CompactTranslatesTensorProductDict, BSplineTensorProductDict, evaluation_operator
export bspline_platform, primal, dual, sampler, A, Zt, ×


export ⊗

include("basisfunctions.jl")

include("operators/tensor_circulant_operator.jl")
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
