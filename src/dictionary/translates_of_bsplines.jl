# translates_of_bsplines.jl
using CardinalBSplines: evaluate_periodic_centered_BSpline, evaluate_periodic_centered_BSpline_derivative
abstract type DiffPeriodicBSplineBasis{T<:Real,K,D} <: PeriodicEquispacedTranslates{T,T}
end

# struct BSplineBasis{T<:Real,K,D}

strings(d::DiffPeriodicBSplineBasis) = (string(d),
    ("length = $(length(d))",
     "$(domaintype(d)) -> $(codomaintype(d))",
     "support = $(support(d))",
     "degree = $(degree(d))"),
     )

# compatible_interpolationgrid(dict::DiffPeriodicBSplineBasis, grid::AbstractEquispacedGrid) =
#     (size(grid) == size(dict)) && (support(grid)≈support(dict)) && (iseven(degree(dict)) ?
#         grid isa MidpointEquispacedGrid :
#         grid isa PeriodicEquispacedGrid)

support(dict::DiffPeriodicBSplineBasis) = UnitInterval{domaintype(dict)}()
measure(dict::DiffPeriodicBSplineBasis{T}) where T = FourierMeasure{T}()
translationgrid(dict::DiffPeriodicBSplineBasis{T}) where T = PeriodicEquispacedGrid(length(dict), zero(T),one(T))

# interpolation_grid(dict::DiffPeriodicBSplineBasis) = isodd(degree(dict)) ?
#     PeriodicEquispacedGrid(length(dict), support(dict)) :
#     MidpointEquispacedGrid(length(dict), support(dict))

hasgrid_transform(dict::DiffPeriodicBSplineBasis, _, grid::AbstractEquispacedGrid) =
    size(dict)==size(grid) && compatible_interpolationgrid(typeof(dict), grid)

eval_kernel(dict::DICT where DICT<:DiffPeriodicBSplineBasis{T,K,D}, x) where {T,K,D} =
    (n = length(dict); sqrt(T(n))*n^D*evaluate_periodic_centered_BSpline_derivative(Val{K}(), Val{D}(), n*x, n, T))

degree(dict::DICT where DICT <: DiffPeriodicBSplineBasis{T,K} ) where {T,K} = K

Bdiff(dict::DICT where DICT <: DiffPeriodicBSplineBasis{T,K,D} ) where {T,K,D} = D


abstract type PeriodicBSplineBasis{T,K} <: DiffPeriodicBSplineBasis{T,K,0}
end

"""
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct BSplineTranslatesBasis{T,K,SCALED} <: PeriodicBSplineBasis{T,K}
    n               :: Int
end


Base.length(dict::BSplineTranslatesBasis) = dict.n
Base.size(dict::BSplineTranslatesBasis) = (length(dict),)

Base.similar(d::BSplineTranslatesBasis{S,K,SCALED}, ::Type{T}, n::Int) where {S,T,K,SCALED} = BSplineTranslatesBasis{T,K,SCALED}(n)

BSplineTranslatesBasis(n::Int, degree::Int, ::Type{T} = Float64; options...) where {T} =
    BSplineTranslatesBasis{T}(n, degree; options...)

BSplineTranslatesBasis{T}(n::Int, degree::Int; options...) where {T} = BSplineTranslatesBasis{T,degree}(n; options...)
BSplineTranslatesBasis{T,degree}(n::Int; scaled=true) where {T,degree} = BSplineTranslatesBasis{T,degree,scaled}(n)



eval_kernel(dict::BSplineTranslatesBasis{T,K,true}, x) where {K,T} =
    (n = length(dict); sqrt(T(n))*evaluate_periodic_centered_BSpline(Val{K}(), n*x, n, T))
eval_kernel(dict::BSplineTranslatesBasis{T,K,false}, x) where {K,T} =
    (n = length(dict); evaluate_periodic_centered_BSpline(Val{K}(), n*x, n, T))

scaled(::BSplineTranslatesBasis{T,K,SCALED}) where {T,K,SCALED} = SCALED

name(dict::BSplineTranslatesBasis) = "Periodic equispaced translates of B spline of degree $(degree(dict))"

==(dict1::BSplineTranslatesBasis{T1,K1}, dict2::BSplineTranslatesBasis{T2,K2}) where {K1,K2,T1,T2} = T1==T2 && K1==K2 && length(dict1)==length(dict2)

instantiate(::Type{BSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = BSplineTranslatesBasis{T}(n,3)

resize(dict::BSplineTranslatesBasis{T}, n::Int) where {T} = BSplineTranslatesBasis{T,degree(dict)}(n)

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
    r = CardinalBSplines.shifted_spline_integral(K, abs(i - 1))
    r += CardinalBSplines.shifted_spline_integral(K, length(dict)-abs(i - 1))
    if SCALED
        convert(T,r)
    else
        convert(T,r)/convert(T,length(dict))
    end
end

"""
  Basis consisting of differentiated dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct DiffBSplineTranslatesBasis{T,K,D} <: DiffPeriodicBSplineBasis{T,K,D}
    n               :: Int
end

Base.length(dict::DiffBSplineTranslatesBasis) = dict.n
Base.size(dict::DiffBSplineTranslatesBasis) = (length(dict),)

Base.similar(d::DiffBSplineTranslatesBasis{S,K,D}, ::Type{T}, n::Int) where {S,T,K,D} =
    DiffBSplineTranslatesBasis{T,K,D}(n)



DiffBSplineTranslatesBasis(n::Int, degree::Int, diff::Int, ::Type{T} = Float64; options...) where {T} =
    DiffBSplineTranslatesBasis{T}(n, degree, diff; options...)

DiffBSplineTranslatesBasis{T}(n::Int, degree::Int, diff::Int; options...) where {T} = DiffBSplineTranslatesBasis{T,degree,diff}(n; options...)
DiffBSplineTranslatesBasis{T,degree}(n::Int; D=0) where {T,degree} = DiffBSplineTranslatesBasis{T,degree,D}(n)

instantiate(::Type{DiffBSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = DiffBSplineTranslatesBasis(n,3,1,T)

resize(dict::DiffBSplineTranslatesBasis, n::Int) = DiffBSplineTranslatesBasis(n, degree(dict), Bdiff(dict), codomaintype(dict))
