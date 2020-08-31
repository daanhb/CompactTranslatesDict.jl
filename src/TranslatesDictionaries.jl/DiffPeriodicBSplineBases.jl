module DiffPeriodicBSplineBases

using ..TranslatesDictionaries, BasisFunctions, DomainSets, GridArrays, ..PeriodicEquispacedTranslatesDicts

using CardinalBSplines: evaluate_centered_BSpline, evaluate_centered_BSpline_derivative,
    evaluate_centered_gauss_BSpline, shifted_spline_integral, evaluate_BSpline, evaluate_BSpline_derivative
using BasisFunctions: difforder, VerticalBandedMatrix

import BasisFunctions: strings, support, measure, size, hasgrid_transform, length, similar,
    resize, innerproduct_native, name, diff, derivative_dict, differentiation, hasderivative
import ..TranslatesDictionaries: translationgrid, unsafe_eval_kernel, kernel_support, eval_kernel, unsafe_eval_kernel_derivative
import ..PeriodicEquispacedTranslatesDicts: firstgramcolumn
import Base: ==


abstract type AbstractDiffPeriodicBSplineBasis{T<:Real,K,D,SCALED,CENTERED} <: PeriodicEquispacedTranslates{T,T,:sum}
end

strings(d::AbstractDiffPeriodicBSplineBasis) = (string(d),
    ("length = $(length(d))",
     "$(domaintype(d)) -> $(codomaintype(d))",
     "support = $(support(d))",
     "degree = $(degree(d)) (diff=$(Bdiff(d)))"),
     )

support(dict::AbstractDiffPeriodicBSplineBasis) =
    UnitInterval{domaintype(dict)}()
measure(dict::AbstractDiffPeriodicBSplineBasis{T}) where T =
    FourierMeasure{T}()
translationgrid(dict::AbstractDiffPeriodicBSplineBasis{T}) where T =
    PeriodicEquispacedGrid(length(dict), zero(T),one(T))
hasgrid_transform(dict::AbstractDiffPeriodicBSplineBasis, _, grid::AbstractEquispacedGrid) =
    size(dict)==size(grid) && compatible_interpolationgrid(dict, grid)

export degree
degree(dict::DICT where DICT <: AbstractDiffPeriodicBSplineBasis{T,K} ) where {T,K} = K
export Bdiff
Bdiff(dict::DICT where DICT <: AbstractDiffPeriodicBSplineBasis{T,K,D} ) where {T,K,D} = D
length(dict::AbstractDiffPeriodicBSplineBasis) = dict.n
size(dict::AbstractDiffPeriodicBSplineBasis) = (length(dict),)
scaled(::AbstractDiffPeriodicBSplineBasis{T,K,D,SCALED}) where {T,K,D,SCALED} = SCALED
centered(::AbstractDiffPeriodicBSplineBasis{T,K,D,S,C}) where {T,K,D,S,C} = C

eval_kernel(dict::AbstractDiffPeriodicBSplineBasis, x) =
    unsafe_eval_kernel(dict, x)
unsafe_eval_kernel(dict::AbstractDiffPeriodicBSplineBasis{T,K,D,true}, x) where {T,K,D} =
    sqrt(T(length(dict)))*unscaled_unsafe_eval_kernel(dict, x)
unsafe_eval_kernel_derivative(dict::AbstractDiffPeriodicBSplineBasis{T,K,D,true}, x, order) where {T,K,D} =
    sqrt(T(length(dict)))*unscaled_unsafe_eval_kernel_derivative(dict, x, order)
unsafe_eval_kernel(dict::AbstractDiffPeriodicBSplineBasis{T,K,D,false}, x)  where {T,K,D} =
    unscaled_unsafe_eval_kernel(dict, x)
unsafe_eval_kernel_derivative(dict::AbstractDiffPeriodicBSplineBasis{T,K,D,false}, x, order) where {T,K,D} =
    unscaled_unsafe_eval_kernel_derivative(dict, x, order)

"""
    struct DiffBSplineTranslatesBasis{T,K,D} <: AbstractDiffPeriodicBSplineBasis{T,K,D}

Basis consisting of differentiated dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct DiffPeriodicBSplineBasis{T,K,D,SCALED,CENTERED} <: AbstractDiffPeriodicBSplineBasis{T,K,D,SCALED,CENTERED}
    n   ::  Int
end
name(dict::DiffPeriodicBSplineBasis) = "Periodic equispaced translates of $(centered(dict) ? "" : "centered ")B spline of degree $(degree(dict))"
similar(d::DiffPeriodicBSplineBasis{S,K,D,SCALED,CENTERED}, ::Type{T}, n::Int) where {S,T,K,D,SCALED,CENTERED} =
    DiffPeriodicBSplineBasis{T,K,D,SCALED,CENTERED}(n)
DiffPeriodicBSplineBasis(n::Int, degree::Int, diff::Int=0, centered::Bool=true, ::Type{T} = Float64; options...) where {T} =
    DiffPeriodicBSplineBasis{T}(n, degree, diff, centered; options...)

DiffPeriodicBSplineBasis{T}(n::Int, degree::Int = 3, diff::Int=0, centered::Bool=true; options...) where {T} =
    DiffPeriodicBSplineBasis{T,degree,diff}(n; options...)
DiffPeriodicBSplineBasis{T,degree,diff}(n::Int; scaled=true, centered=true) where {T,degree,diff} =
    DiffPeriodicBSplineBasis{T,degree,diff,scaled,centered}(n)

unscaled_unsafe_eval_kernel(dict::DiffPeriodicBSplineBasis{T,K,D,S,true}, x) where {T,K,D,S} =
    (n = length(dict); n^D*evaluate_centered_BSpline_derivative(Val{K}(), Val{D}(), n*x, T))
unscaled_unsafe_eval_kernel_derivative(dict::DiffPeriodicBSplineBasis{T,K,D,true}, x, order::Int) where {T,K,D,S} =
    (n = length(dict); n^(D+order)*evaluate_centered_BSpline_derivative(Val{K}(), Val(D+order), n*x, T))

_shift(K,T) = iseven(K) ? zero(T) : one(T)/2
unscaled_unsafe_eval_kernel(dict::DiffPeriodicBSplineBasis{T,K,D,S,false}, x) where {T,K,D,S} =
    (n = length(dict); n^D*evaluate_BSpline_derivative(Val{K}(), Val{D}(), n*x+_shift(K,T), T))
unscaled_unsafe_eval_kernel_derivative(dict::DiffPeriodicBSplineBasis{T,K,D,false}, x, order::Int) where {T,K,D,S} =
    (n = length(dict); n^(D+order)*evaluate_BSpline_derivative(Val{K}(), Val(D+order), n*x+_shift(K,T), T))

kernel_support(dict::DiffPeriodicBSplineBasis{T,K,D,S,true}) where {T,K,D,S} =
    DomainSets.Interval{:closed,:open}(-convert(coefficienttype(dict),degree(dict)+1)/2/length(dict),convert(coefficienttype(dict),degree(dict)+1)/2/length(dict))
kernel_support(dict::DiffPeriodicBSplineBasis{T,K,D,S,false}) where {T,K,D,S} =
    DomainSets.Interval{:closed,:open}(-_shift(degree(dict),T)/length(dict),convert(coefficienttype(dict),degree(dict)+1-_shift(degree(dict),T))/length(dict))

diff(dict::DiffPeriodicBSplineBasis{T,K,D,SCALED,C}, order; options...) where {T,K,D,SCALED,C} =
    DiffPeriodicBSplineBasis{T,K,D+order,SCALED,C}(length(dict))

hasderivative(::DiffPeriodicBSplineBasis) = true
function derivative_dict(dict::DiffPeriodicBSplineBasis{T,K,D,SCALED,C}, order; options...) where {T,K,D,SCALED,C}
    @assert order+Bdiff(dict) <= degree(dict)
    if isodd(order+Bdiff(dict))
        # knots of derivative_dict do not align with dict
        return DiffPeriodicBSplineBasis{T,K-D-order,0,SCALED,!(C)}(length(dict))
    end
    DiffPeriodicBSplineBasis{T,K-D-order,0,SCALED,C}(length(dict))
end

function differentiation(::Type{T}, src::DiffPeriodicBSplineBasis, dest::DiffPeriodicBSplineBasis, order::Int; options...) where {T}
    @assert Bdiff(dest)==0
    @assert degree(dest)==degree(src)-Bdiff(src)-order
    isodd(order+Bdiff(src)) && (centered(dest) == centered(src)) &&
        error("Alignment of both dictionaries lead to inconsistency")
    length(dest)!=length(src) &&
        error("Dictionaries must have same length")


    N = length(src)
    s = (scaled(src) == scaled(dest)) ?
        one(T) :
        sqrt(T(N))
    scaled(dest) && (s = one(T)/s)

    ArrayOperator(VerticalBandedMatrix(N,N,s*N^order*[(-1)^k*binomial(order, k) for k in 0:order],1,-order),
        src, dest)
end

==(::DiffPeriodicBSplineBasis, ::DiffPeriodicBSplineBasis) = false
==(dict1::DiffPeriodicBSplineBasis{T,K,D,S,C}, dict2::DiffPeriodicBSplineBasis{T,K,D,S,C}) where {T,K,S,D,C} = length(dict1)==length(dict2)

resize(dict::DiffPeriodicBSplineBasis{T}, n::Int) where {T} =
    DiffPeriodicBSplineBasis{T,degree(dict),Bdiff(dict),scaled(dict),centered(dict)}(n)

export BSplineTranslatesBasis
"""
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
const BSplineTranslatesBasis{T,K,SCALED,CENTERED}  = DiffPeriodicBSplineBasis{T,K,0,SCALED,CENTERED}

similar(d::BSplineTranslatesBasis{S,K,SCALED,CENTERED}, ::Type{T}, n::Int) where {S,T,K,SCALED,CENTERED} =
    BSplineTranslatesBasis{T,K,SCALED,CENTERED}(n)
BSplineTranslatesBasis(n::Int, degree::Int, ::Type{T} = Float64; options...) where {T} =
    BSplineTranslatesBasis{T}(n, degree; options...)

BSplineTranslatesBasis{T}(n::Int, degree::Int = 3; options...) where {T} =
    BSplineTranslatesBasis{T,degree}(n; options...)
BSplineTranslatesBasis{T,degree}(n::Int; scaled=true,centered=true) where {T,degree} =
    DiffPeriodicBSplineBasis{T,degree,0,scaled,centered}(n)

function BSplineTranslatesBasis(n::Int, degree::Int, a::T, b::T, ::Type{S}=T; options...) where {T,S}
    Q = promote_type(T, S, typeof(a/b))
    rescale(BSplineTranslatesBasis{Q}(n, degree;options...), Q(a), Q(b))
end

unscaled_unsafe_eval_kernel(dict::BSplineTranslatesBasis{T,K,S,true}, x) where {K,T,S} =
    (n = length(dict); evaluate_centered_BSpline(Val{K}(), n*x, T))
unscaled_unsafe_eval_kernel(dict::BSplineTranslatesBasis{T,K,S,false}, x) where {K,T,S} =
    (n = length(dict); evaluate_BSpline(Val{K}(), n*x+_shift(K,T), T))

innerproduct_native(d1::BSplineTranslatesBasis, i, d2::BSplineTranslatesBasis, j, measure::FourierMeasure; warnslow=false, options...)  =
     BasisFunctions.default_dict_innerproduct(d1, i, d2, j, measure; warnslow=warnslow, options...)

function firstgramcolumn(dict::BSplineTranslatesBasis, measure::FourierMeasure; T=coefficienttype(s, domaintype(measure)), options...)
    firstcolumn = zeros(T, length(dict))
    for i in 1:length(dict)
        firstcolumn[i] = firstgramcolumnelement(dict, measure, i; T=T, options...)
    end
    firstcolumn
end


@inline function firstgramcolumnelement(dict::BSplineTranslatesBasis{S,K,SCALED}, measure::FourierMeasure, i::Int; T=coefficienttype(s), options...) where {S,K,SCALED}
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


"""
  Basis consisting of dilated, translated, and periodized gauss clocks similar to cardinal B splines on the interval [0,1].
"""
struct GaussTranslatesBasis{T,K,SCALED} <: AbstractDiffPeriodicBSplineBasis{T,K,0,SCALED,true}
    n   :: Int
end
name(dict::GaussTranslatesBasis) = "Periodic equispaced translates of gaussians similar to B spline of degree $(degree(dict))"
similar(d::GaussTranslatesBasis{S,K,SCALED}, ::Type{T}, n::Int) where {S,T,K,SCALED} = GaussTranslatesBasis{T,K,SCALED}(n)
GaussTranslatesBasis(n::Int, degree::Int, ::Type{T} = Float64; options...) where {T} =
    GaussTranslatesBasis{T}(n, degree; options...)

GaussTranslatesBasis{T}(n::Int, degree::Int = 3; options...) where {T} = GaussTranslatesBasis{T,degree}(n; options...)
GaussTranslatesBasis{T,degree}(n::Int; scaled=true) where {T,degree} = GaussTranslatesBasis{T,degree,scaled}(n)

function GaussTranslatesBasis(n::Int, degree::Int, a::T, b::T, ::Type{S}=T; options...) where {T,S}
    Q = promote_type(T, S, typeof(a/b))
    rescale(GaussTranslatesBasis{Q}(n, degree;options...), Q(a), Q(b))
end

unscaled_unsafe_eval_kernel(dict::GaussTranslatesBasis{T,K}, x) where {K,T} =
    (n = length(dict); evaluate_centered_gauss_BSpline(Val{K}(), n*x, T))


==(::GaussTranslatesBasis, ::GaussTranslatesBasis) = false
==(dict1::GaussTranslatesBasis{T,K,S}, dict2::GaussTranslatesBasis{T,K,S}) where {T,S,K} = length(dict1)==length(dict2)

resize(dict::GaussTranslatesBasis{T}, n::Int) where {T} = GaussTranslatesBasis{T,degree(dict),scaled(dict)}(n)

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


end
