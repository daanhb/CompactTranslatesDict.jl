module CompactPeriodicEquispacedTranslatesDuals

using BasisFunctions, DomainSets, ..TranslatesDictionaries, ..PeriodicEquispacedTranslatesDicts

using BasisFunctions: VerticalBandedMatrix
using InfiniteVectors: downsample, subvector, PeriodicInfiniteVector, sublength, shift!
using ..CompactInfiniteVectors: compactinfinitevector
using GridArrays: PeriodicEquispacedGrid

import BasisFunctions: string, strings, name, unsafe_eval_element, grid_evaluation_operator, size, length, support
import CompactTranslatesDict: eval_kernel
import CardinalBSplines: minimalK

export CompactPeriodicEquispacedTranslatesDual
"""
    struct CompactPeriodicEquispacedTranslatesDual{T,K,DICT<:PeriodicEquispacedTranslates{T,K,PERIODIZATION}}
            <: PeriodicEquispacedTranslates{T,K,PERIODIZATION}

Compact, discrete, periodic dual --- with respect to a `DiracCombMeasure` ---
to the a dictionary consisting of equispaced
periodic translates of a kernel (see `PeriodicEquispacedTranslates` of `CompactTranslatesDict`).

Since this basis is discrete, it can only be evaluated in an `PeriodicEquispacedGrid`.

# Example
```jldocs
julia> using BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals, CompactTranslatesDict

julia> B = CompactPeriodicEquispacedTranslatesDual(BSplineTranslatesBasis(10,3,-1,1), 2)
Equispaced translates of a discrete kernel dual
    ↳ Periodic equispaced translates of a periodic kernel function
      ↳ length = 10
      ↳ Float64 -> Float64
      ↳ support = -1.0..1.0
    ↳ m = 2

julia> A = evaluation_operator(B, PeriodicEquispacedGrid(10, BasisFunctions.support(B)))
Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}
```
"""
struct CompactPeriodicEquispacedTranslatesDual{T,S,PERIODIZATION,DICT<:Dictionary1d} <: PeriodicEquispacedTranslates{T,S,PERIODIZATION}
    dict   ::   DICT
    m      ::   Int
    minimalK::Int
end

function CompactPeriodicEquispacedTranslatesDual(dict::Dictionary1d{T,S}, m::Int) where {T,S}
    N =  m*length(dict)
    b = signal(dict, m)
    minK = -1
    for K in 0:2sublength(b)
        try
            inv(b,m,K=K)
            minK = K
            break;
        catch e
            if e isa ErrorException && e.msg[1:20] == "Can not find compact"
                nothing
            else
                throw(e)
            end
        end
    end
    if minK==-1
        error("No compact dual found")
    end
    CompactPeriodicEquispacedTranslatesDual{T,S,:sum,typeof(dict)}(dict, m, minK)
end
size(dict::CompactPeriodicEquispacedTranslatesDual) = (length(dict.dict),)
length(dict::CompactPeriodicEquispacedTranslatesDual) = length(dict.dict)
support(dict::CompactPeriodicEquispacedTranslatesDual) = support(dict.dict)

name(dict::CompactPeriodicEquispacedTranslatesDual)= "Equispaced translates of a discrete kernel dual"

strings(d::CompactPeriodicEquispacedTranslatesDual) = (string(d),strings(d.dict),(
     "m = $(d.m)",),
     )

unsafe_eval_element(dict::CompactPeriodicEquispacedTranslatesDual, i, x) =
    error("`CompactPeriodicEquispacedTranslatesDual` can only be evaluated in `PeriodicEquispacedGrid`")
grid_evaluation_operator(dict::CompactPeriodicEquispacedTranslatesDual, gb::GridBasis, grid::AbstractGrid) =
    error("`CompactPeriodicEquispacedTranslatesDual` can only be evaluated in `PeriodicEquispacedGrid`")
grid_evaluation_operator(dict::CompactPeriodicEquispacedTranslatesDual, dgs::GridBasis, grid::AbstractEquispacedGrid) =
    error("`CompactPeriodicEquispacedTranslatesDual` can only be evaluated in `PeriodicEquispacedGrid`")
eval_kernel(dict::CompactPeriodicEquispacedTranslatesDual, x) =
    error("`CompactPeriodicEquispacedTranslatesDual` can only be evaluated in `PeriodicEquispacedGrid`")

function grid_evaluation_operator(dict::CompactPeriodicEquispacedTranslatesDual{T}, gb::GridBasis, grid::PeriodicEquispacedGrid;
        threshold=100eps(T), options...) where {T}
    @assert support(dict) ≈ support(grid)
    m_dict = dict.m
    m_frac = m_dict*length(dict) / length(grid)

    @assert round(m_frac) ≈ m_frac
    m_frac = round(Int, m_frac)

    N =  m_dict*length(dict)
    b = signal(dict.dict, m_dict)
    dual_signal = PeriodicInfiniteVector(inv(b, m_dict, threshold; K=dict.minimalK)', N)

    v = downsample(dual_signal, m_frac)

    # Convert PeriodicInfiniteVector to VerticalBandedMatrix
    mask = 0 .==subvector(v)
    a, b = findfirst(mask), findlast(mask)

    if a!=nothing && b!= nothing
        if a == 1 || b==length(mask)
            bw = findlast(.!(mask)), findfirst(.!(mask))
            a = subvector(v)[bw[2]:bw[1]]
        else
            bw = a-1, b+1
            a = [subvector(v)[bw[2]:end];subvector(v)[1:bw[1]]]
        end
    else
        bw = 0, 1
        a = subvector(v)
    end

    (length(grid)<length(dict)) && (error("Not implemented."))
    shift = length(grid) / length(dict)
    (round(shift) ≈ shift) || error("not implemented")
    shift = round(Int, shift)

    ArrayOperator(VerticalBandedMatrix(length(grid), length(dict), a, shift, -1+bw[2]-length(grid)), dict, GridBasis{coefficienttype(dict)}(grid))
end

signal(dict::Dictionary, m::Int) =
    shift!(compactinfinitevector(dict, PeriodicEquispacedGrid(m*length(dict),support(dict)),),-1)

end
