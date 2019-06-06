__precompile__()
module CompactTranslatesDict
using Reexport, CardinalBSplines, DomainSets, FillArrays, LinearAlgebra, FFTW
@reexport using BasisFunctions

import Base: ==, getindex, length, size, unsafe_getindex, checkbounds, step

# Dictionaries
using BasisFunctions
using BasisFunctions: ShiftedIndex, ShiftedIndexList, op_eltype, dense_evaluation_operator, Measure
using BasisFunctions: UniformDiracCombMeasure, grid, default_gramoperator
import BasisFunctions: hasmeasure, measure, gramoperator, innerproduct_native, mixedgramoperator
import BasisFunctions: support, interpolation_grid, period
import BasisFunctions: approx_length
import BasisFunctions: transform_from_grid, transform_to_grid, hasgrid_transform, iscompatible
import BasisFunctions: grid_evaluation_operator, unsafe_eval_element, similar_dictionary
import BasisFunctions: matrix, apply!, adjoint, unsafe_wrap_operator, matrix!
import BasisFunctions: extension_operator, restriction_operator, resize
import BasisFunctions: elements, native_index, src, dest, *, +, -, inv, pinv, element
import BasisFunctions: instantiate, isbasis, isbiorthogonal, name, ordering



export BSplineTranslatesBasis, DiffBSplineTranslatesBasis
export degree, Bdiff, src, dest
export CompactTranslatesTensorProductDict, BSplineTensorProductDict, evaluation_operator



# include("operators/tensor_circulant_operator.jl")
include("dictionary/translates.jl")
include("dictionary/translates_of_bsplines.jl")
include("dictionary/tensor.jl")
end
