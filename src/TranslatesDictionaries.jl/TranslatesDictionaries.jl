module TranslatesDictionaries

using BasisFunctions, DomainSets, GridArrays, Reexport

using BasisFunctions: DomainLebesgueMeasure

import BasisFunctions: support, size, length, ordering, unsafe_eval_element, hasinterpolationgrid,
    interpolation_grid, hasmeasure, measure, name,
    isperiodic, period, similar, rescale, dict_in_support
import Base: step

export Translates
"""
    abstract type Translates{T,S} <: Dictionary{T,S}
"""
abstract type Translates{T,S} <: Dictionary{T,S}
end

compatible_translationgrid(dict::Translates,grid::AbstractGrid) =
    size(dict) == size(grid) && compatible_translationgrid(typeof(dict), grid)
compatible_translationgrid(::Type{<:Translates}, grid::AbstractGrid) =
    first(grid) ∈ coverdomain(grid)

export translationgrid
"""
    translationgrid(dict::Translates)

The translation grid of a Translates `dict`.
"""
translationgrid(dict::Translates) = dict.grid
support(dict::Translates) = coverdomain(translationgrid(dict))
support(dict::Translates, idx) = support(dict)
size(dict::Translates) = size(translationgrid(dict))
length(dict::Translates) = length(translationgrid(dict))
ordering(dict::Translates) = eachindex(translationgrid(dict))
unsafe_eval_element(dict::Translates, idxn, x) =
    unsafe_eval_kernel(dict, x-translationgrid(dict)[idxn])
dict_in_support(dict::Translates, idxn, x) =
    (x-translationgrid(dict)[idxn]) ∈ kernel_support(dict)

export kernel_support
"""
    kernel_support(dict::Translates)

The support of the kernel of a Translates dictionary
"""
kernel_support(dict::Translates) = dict.kernel_support
hasinterpolationgrid(::Translates) = true
interpolation_grid(dict::Translates) = translationgrid(dict)
compatible_interpolationgrid(dict::Translates, grid::AbstractGrid) =
    support(dict)≈coverdomain(grid) && compatible_translationgrid(dict, grid)
hasmeasure(dict::Translates) = true
measure(dict::Translates) = DomainLebesgueMeasure(support(dict))


export eval_kernel
"""
    eval_kernel(dict::Translates, x)

Evaluate the kernel of `dict`in `x`.
"""
function eval_kernel(dict::Translates{T}, x) where {T}
    if x ∈ kernel_support(dict)
        dict.kernel(x)
    else
        zero(T)
    end
end

unsafe_eval_kernel(dict::Translates, x) =
    dict.kernel(x)

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
struct GenericTranslates{T,S,GRID<:AbstractGrid,DOMAIN<:Domain,KERN} <: Translates{T,S}
    grid::GRID
    kernel::KERN
    kernel_support::DOMAIN
    function GenericTranslates(grid, kernel, kernel_support=DomainSets.FullSpace{typeof(first(grid))}())
        T = eltype(grid)
        S = typeof(kernel(first(grid)))
        @assert compatible_translationgrid(GenericTranslates, grid)
        new{T,S,typeof(grid),typeof(kernel_support),typeof(kernel)}(grid, kernel, kernel_support)
    end
end

export EquispacedTranslates
"""
    abstract type EquispacedTranslates{T,S} <: Translates{T,S}
"""
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
struct GenericEquispacedTranslates{T,S,GRID<:AbstractGrid,DOMAIN<:Domain,KERN} <: EquispacedTranslates{T,S}
    grid::GRID
    kernel::KERN
    kernel_support::DOMAIN
    function GenericEquispacedTranslates(grid, kernel, kernel_support=DomainSets.FullSpace{typeof(first(grid))}())
        T = eltype(grid)
        S = typeof(kernel(first(grid)))
        @assert compatible_translationgrid(GenericEquispacedTranslates, grid)
        new{T,S,typeof(grid),typeof(kernel_support),typeof(kernel)}(grid, kernel, kernel_support)
    end
end

end
