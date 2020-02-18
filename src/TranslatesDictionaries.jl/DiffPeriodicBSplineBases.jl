module DiffPeriodicBSplineBases

using ..TranslatesDictionaries, BasisFunctions, DomainSets, GridArrays, ..PeriodicEquispacedTranslatesDicts

using CardinalBSplines: evaluate_centered_BSpline, evaluate_centered_BSpline_derivative,
    evaluate_centered_gauss_BSpline, shifted_spline_integral

import BasisFunctions: strings, support, measure, size, hasgrid_transform, length, similar,
    instantiate, resize, innerproduct_native, name
import ..TranslatesDictionaries: translationgrid, unsafe_eval_kernel, kernel_support, eval_kernel
import ..PeriodicEquispacedTranslatesDicts: firstgramcolumn
import Base: ==


abstract type DiffPeriodicBSplineBasis{T<:Real,K,D} <: PeriodicEquispacedTranslates{T,T,:sum}
end

strings(d::DiffPeriodicBSplineBasis) = (string(d),
    ("length = $(length(d))",
     "$(domaintype(d)) -> $(codomaintype(d))",
     "support = $(support(d))",
     "degree = $(degree(d))"),
     )

support(dict::DiffPeriodicBSplineBasis) =
    UnitInterval{domaintype(dict)}()
measure(dict::DiffPeriodicBSplineBasis{T}) where T =
    FourierMeasure{T}()
translationgrid(dict::DiffPeriodicBSplineBasis{T}) where T =
    PeriodicEquispacedGrid(length(dict), zero(T),one(T))
hasgrid_transform(dict::DiffPeriodicBSplineBasis, _, grid::AbstractEquispacedGrid) =
    size(dict)==size(grid) && compatible_interpolationgrid(dict, grid)
eval_kernel(dict::DICT where DICT<:DiffPeriodicBSplineBasis{T,K,D}, x) where {T,K,D} =
    unsafe_eval_kernel(dict, x)
unsafe_eval_kernel(dict::DICT where DICT<:DiffPeriodicBSplineBasis{T,K,D}, x) where {T,K,D} =
    (n = length(dict); sqrt(T(n))*n^D*evaluate_centered_BSpline_derivative(Val{K}(), Val{D}(), n*x, T))
kernel_support(dict::DiffPeriodicBSplineBasis) =
    DomainSets.Interval{:closed,:open}(-convert(coefficienttype(dict),degree(dict)+1)/2/length(dict),convert(coefficienttype(dict),degree(dict)+1)/2/length(dict))

export degree
degree(dict::DICT where DICT <: DiffPeriodicBSplineBasis{T,K} ) where {T,K} = K
export Bdiff
Bdiff(dict::DICT where DICT <: DiffPeriodicBSplineBasis{T,K,D} ) where {T,K,D} = D

export PeriodicBSplineBasis
"""
    abstract type PeriodicBSplineBasis{T,K} <: DiffPeriodicBSplineBasis{T,K,0}
"""
abstract type PeriodicBSplineBasis{T,K} <: DiffPeriodicBSplineBasis{T,K,0}
end

length(dict::PeriodicBSplineBasis) = dict.n
size(dict::PeriodicBSplineBasis) = (length(dict),)

for (TYPE,fun) in zip(
        (:BSplineTranslatesBasis, :GaussTranslatesBasis),
        (:evaluate_centered_BSpline, :evaluate_centered_gauss_BSpline)
    )
    @eval begin

        struct $(TYPE){T,K,SCALED} <: PeriodicBSplineBasis{T,K}
            n               :: Int
        end

        similar(d::$(TYPE){S,K,SCALED}, ::Type{T}, n::Int) where {S,T,K,SCALED} = $(TYPE){T,K,SCALED}(n)
        $(TYPE)(n::Int, degree::Int, ::Type{T} = Float64; options...) where {T} =
            $(TYPE){T}(n, degree; options...)

        $(TYPE){T}(n::Int, degree::Int; options...) where {T} = $(TYPE){T,degree}(n; options...)
        $(TYPE){T,degree}(n::Int; scaled=true) where {T,degree} = $(TYPE){T,degree,scaled}(n)

        function $(TYPE)(n::Int, degree::Int, a::T, b::T, ::Type{S}=T; options...) where {T,S}
            Q = promote_type(T, S, typeof(a/b))
            rescale($(TYPE){Q}(n, degree;options...), Q(a), Q(b))
        end

        unsafe_eval_kernel(dict::$(TYPE){T,K,true}, x) where {K,T} =
            (n = length(dict); sqrt(T(n))*$(fun)(Val{K}(), n*x, T))
        unsafe_eval_kernel(dict::$(TYPE){T,K,false}, x) where {K,T} =
            (n = length(dict); $(fun)(Val{K}(), n*x, T))

        scaled(::$(TYPE){T,K,SCALED}) where {T,K,SCALED} = SCALED

        ==(dict1::$(TYPE){T1,K1}, dict2::$(TYPE){T2,K2}) where {K1,K2,T1,T2} = T1==T2 && K1==K2 && length(dict1)==length(dict2)

        instantiate(::Type{$(TYPE)}, n::Int, ::Type{T}) where {T} = $(TYPE){T}(n,3)

        resize(dict::$(TYPE){T}, n::Int) where {T} = $(TYPE){T,degree(dict)}(n)

        innerproduct_native(d1::$(TYPE), i, d2::$(TYPE), j, measure::FourierMeasure; warnslow=false, options...)  =
             BasisFunctions.default_dict_innerproduct(d1, i, d2, j, measure; warnslow=warnslow, options...)

        function firstgramcolumn(dict::$(TYPE), measure::FourierMeasure; T=coefficienttype(s, domaintype(measure)), options...)
            firstcolumn = zeros(T, length(dict))
            for i in 1:length(dict)
                firstcolumn[i] = firstgramcolumnelement(dict, measure, i; T=T, options...)
            end
            firstcolumn
        end


        @inline function firstgramcolumnelement(dict::$(TYPE){S,K,SCALED}, measure::FourierMeasure, i::Int; T=coefficienttype(s), options...) where {S,K,SCALED}
            if length(dict) <= 2K+1
                return convert(T, innerproduct(dict, ordering(dict)[i], dict, ordering(dict)[1], measure; options...))
            end
            r = shifted_spline_integral(K, abs(i - 1))
            r += shifted_spline_integral(K, length(dict)-abs(i - 1))
            if SCALED
                convert(T,r)
            else
                convert(T,r)/convert(T,length(dict))
            end
        end
    end
end

export BSplineTranslatesBasis, GaussTranslatesBasis
@doc """
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
BSplineTranslatesBasis

@doc """
  Basis consisting of dilated, translated, and periodized gauss clocks similar to cardinal B splines on the interval [0,1].
"""
GaussTranslatesBasis

name(dict::BSplineTranslatesBasis) = "Periodic equispaced translates of B spline of degree $(degree(dict))"
name(dict::GaussTranslatesBasis) = "Periodic equispaced translates of gaussians similar to B spline of degree $(degree(dict))"

kernel_support(dict::GaussTranslatesBasis{T}, threshold = eps(T)) where T =
    approximate_support(dict, threshold)

export approximate_kernel_support
"""
    approximate_kernel_support(dict::Translates, threshold = eps(T))

The support where a kernel is larger than a threshold
"""
function approximate_kernel_support(dict::GaussTranslatesBasis{T}, threshold = eps(T)) where {T}
    a = scaled(dict) ?
        sqrt(-log(threshold*sqrt(T(degree(dict))/6/length(dict)))*T(degree(dict)+1)/6 )/length(dict) :
        sqrt(-log(threshold*sqrt(T(degree(dict))/6             ))*T(degree(dict)+1)/6 )/length(dict)
    -a..a
end

approximate_kernel_support(dict::Translates, threshold = eps(T)) where {T} =
    kernel_support(dict)

export DiffBSplineTranslatesBasis
"""
    struct DiffBSplineTranslatesBasis{T,K,D} <: DiffPeriodicBSplineBasis{T,K,D}

Basis consisting of differentiated dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct DiffBSplineTranslatesBasis{T,K,D} <: DiffPeriodicBSplineBasis{T,K,D}
    n               :: Int
end

length(dict::DiffBSplineTranslatesBasis) = dict.n
size(dict::DiffBSplineTranslatesBasis) = (length(dict),)

similar(d::DiffBSplineTranslatesBasis{S,K,D}, ::Type{T}, n::Int) where {S,T,K,D} =
    DiffBSplineTranslatesBasis{T,K,D}(n)



DiffBSplineTranslatesBasis(n::Int, degree::Int, diff::Int, ::Type{T} = Float64; options...) where {T} =
    DiffBSplineTranslatesBasis{T}(n, degree, diff; options...)

DiffBSplineTranslatesBasis{T}(n::Int, degree::Int, diff::Int; options...) where {T} = DiffBSplineTranslatesBasis{T,degree,diff}(n; options...)
DiffBSplineTranslatesBasis{T,degree}(n::Int; D=0) where {T,degree} = DiffBSplineTranslatesBasis{T,degree,D}(n)

instantiate(::Type{DiffBSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = DiffBSplineTranslatesBasis(n,3,1,T)

resize(dict::DiffBSplineTranslatesBasis, n::Int) = DiffBSplineTranslatesBasis(n, degree(dict), Bdiff(dict), codomaintype(dict))


export CompactTranslatesTensorProductDict, BSplineTensorProductDict
const CompactTranslatesTensorProductDict{N,TUPLE,S,T} = BasisFunctions.TensorProductDict{N,TUPLE,S,T} where {N,TUPLE<: Tuple{Vararg{CompactTranslationDict}},S,T}
const BSplineTensorProductDict{N,TUPLE,S,T} = BasisFunctions.TensorProductDict{N,TUPLE,S,T} where {N,TUPLE<: Tuple{Vararg{DiffPeriodicBSplineBasis}},S,T}

degree(T::BSplineTensorProductDict) = map(degree, elements(T))
Bdiff(T::BSplineTensorProductDict) = map(Bdiff, elements(T))

evaluation(s::CompactTranslatesTensorProductDict, dgs::GridBasis, grid::ProductGrid; options...) =
    tensorproduct([evaluation(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)

tensor_product_dict(size::NTuple{N,Int}, degree::NTuple{N,Int}, diff::NTuple{N,Int}=ntuple(k->0,Val(N)),::Type{T}=Float64) where {N,T} =
    TensorProductDict([DiffBSplineTranslatesBasis(size[i], degree[i], diff[i], T) for i in 1:length(size)]...)

end
