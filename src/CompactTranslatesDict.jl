__precompile__()
module CompactTranslatesDict
using Reexport, CardinalBSplines, DomainSets, FillArrays, LinearAlgebra, FFTW, BasisFunctions

import Base: ==, getindex, length, size, unsafe_getindex, checkbounds, step, similar

# Dictionaries
import BasisFunctions: support, name, string, strings, isperiodic, period, gramoperator,grid_evaluation_operator,
    transform_from_grid, hasgrid_transform, ordering, measure, unsafe_eval_element, hasmeasure, interpolation_grid,
    hasinterpolationgrid, extension_operator, restriction_operator, instantiate, resize, rescale

using BasisFunctions: Measure, GenericLebesgueMeasure, default_mixedgramoperator_discretemeasure,
    DiscreteMeasure, default_gramoperator, op_eltype
using BasisFunctions.GridArrays: similargrid


export BSplineTranslatesBasis, DiffBSplineTranslatesBasis
export degree, Bdiff, src, dest
export CompactTranslatesTensorProductDict, BSplineTensorProductDict, evaluation_operator

include("dictionary/periodic_interval.jl")
abstract type Translates{T,S} <: Dictionary{T,S}
end

compatible_translationgrid(dict::Translates,grid::AbstractGrid) =
    size(dict) == size(grid) && compatible_translationgrid(typeof(dict), grid)
compatible_translationgrid(::Type{<:Translates}, grid::AbstractGrid) =
    first(grid) ∈ support(grid)

export translationgrid, support, eval_kernel
"""
    translationgrid(dict::Translates)

Return the translation grid of a Translates `dict`.
"""
translationgrid(dict::Translates) = dict.grid
support(dict::Translates) = support(translationgrid(dict))
support(dict::Translates, idx) = support(dict, native_index(dict, idx))
size(dict::Translates) = size(translationgrid(dict))
length(dict::Translates) = length(translationgrid(dict))
ordering(dict::Translates) = eachindex(translationgrid(dict))
unsafe_eval_element(dict::Translates, idxn, x) =
    eval_kernel(dict, x-translationgrid(dict)[idxn])
kernel_support(dict::Translates) = dict.kernel_support
hasinterpolationgrid(::Translates) = true
interpolation_grid(dict::Translates) = translationgrid(dict)
compatible_interpolationgrid(dict::Translates, grid::AbstractGrid) =
    support(dict)≈support(grid) && compatible_translationgrid(dict, grid)
hasmeasure(dict::Translates) = true
measure(dict::Translates) = GenericLebesgueMeasure(support(dict))

"""
    eval_kernel(dict::Translates, x)

Evaluate the kernel of `dict`in `x`.
"""
eval_kernel(dict::Translates, x) = dict.kernel(x)
name(::Translates) = "Translates of a kernel function"

export GenericTranslates
"""
struct GenericTranslates{T,S} <: Translates{T,S}

A `Translates` with a general translation grid and kernel

# Example
```jldocs
julia> GenericTranslates(EquispacedGrid(10,0,1), exp)
GenericTranslates
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0
```
"""
struct GenericTranslates{T,S} <: Translates{T,S}
    grid::AbstractGrid
    kernel
    kernel_support::Domain
    function GenericTranslates(grid, kernel, kernel_support=DomainSets.FullSpace{typeof(first(grid))}())
        T = eltype(grid)
        S = typeof(kernel(first(grid)))
        @assert compatible_translationgrid(GenericTranslates, grid)
        new{T,S}(grid, kernel, kernel_support)
    end
end


abstract type EquispacedTranslates{T,S} <: Translates{T,S}
end

compatible_translationgrid(::Type{<:EquispacedTranslates}, grid::AbstractGrid) =
    compatible_translationgrid(Translates, grid) && step(grid) >= 0
step(dict::EquispacedTranslates) = step(translationgrid(dict))
name(::EquispacedTranslates) = "Equispaced translates of a kernel function"

export GenericEquispacedTranslates
"""
struct GenericEquispacedTranslates{T,S} <: Translates{T,S}

A `Translates` with a general equispaced translation grid and kernel

# Example
```jldocs
julia> GenericEquispacedTranslates(EquispacedGrid(10,0,1), exp)
GenericEquispacedTranslates
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0
```
"""
struct GenericEquispacedTranslates{T,S} <: EquispacedTranslates{T,S}
    grid::AbstractGrid
    kernel
    kernel_support::Domain
    function GenericEquispacedTranslates(grid, kernel, kernel_support=DomainSets.FullSpace{typeof(first(grid))}())
        T = eltype(grid)
        S = typeof(kernel(first(grid)))
        @assert compatible_translationgrid(GenericEquispacedTranslates, grid)
        new{T,S}(grid, kernel, kernel_support)
    end
end

"""
    abstract type PeriodicEquispacedTranslates{T,S} <: EquispacedTranslates{T,S}

A dictionary consisting of the equispaced translates of a periodic kernel.
The kernel is periodized during evaluation.
"""
abstract type PeriodicEquispacedTranslates{T,S} <: EquispacedTranslates{T,S}
end

compatible_translationgrid(::Type{<:PeriodicEquispacedTranslates}, grid::AbstractGrid) =
    compatible_translationgrid(EquispacedTranslates, grid) && isperiodic(grid)
isperiodic(::PeriodicEquispacedTranslates) = true
export period
"""
    period(dict::PeriodicEquispacedTranslates)

The period of the periodic dictionary
"""
period(dict::PeriodicEquispacedTranslates) = DomainSets.width(support(dict))
support(dict::PeriodicEquispacedTranslates, idx) =
	PeriodicInterval(translationgrid(dict)[idx]+kernel_support(dict), support(dict))

name(::PeriodicEquispacedTranslates) = "Periodic equispaced translates of a periodic kernel function"
include("periodicequispacedtranslates.jl")


export GenericPeriodicEquispacedTranslates
"""
struct GenericPeriodicEquispacedTranslates{T,S} <: Translates{T,S}

A `Translates` with a general periodic equispaced grid and periodic kernel

# Example
```jldocs
julia> GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(10,0,2π), cos)
GenericPeriodicEquispacedTranslates
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..6.283185307179586
```
"""
struct GenericPeriodicEquispacedTranslates{T,S} <: PeriodicEquispacedTranslates{T,S}
    grid::AbstractGrid
    kernel
    kernel_support::Domain
    function GenericPeriodicEquispacedTranslates(grid, kernel, kernel_support=DomainSets.FullSpace{typeof(first(grid))}())
        T = eltype(grid)
        S = typeof(kernel(first(grid)))
        @assert compatible_translationgrid(GenericPeriodicEquispacedTranslates, grid)
        new{T,S}(grid, kernel, kernel_support)
    end
end

similar(dict::GenericPeriodicEquispacedTranslates{K,S}, ::Type{T}, n) where {K,S,T} =
    GenericPeriodicEquispacedTranslates(similargrid(translationgrid(dict),real(T), n), dict.kernel, dict.kernel_support)

function rescale(dict::PeriodicEquispacedTranslates, a, b)
    map = interval_map(extrema(support(dict))...,   a, b)
    GenericPeriodicEquispacedTranslates(mapped_grid(translationgrid(dict),map),
        x->eval_kernel(dict, inv(map)*x),
        map*kernel_support(dict))
end

include("dictionary/translates_of_bsplines.jl")
include("dictionary/tensor.jl")
end
