# translates_of_bsplines.jl
using CardinalBSplines: evaluate_periodic_Bspline, evaluate_periodic_Bspline_derivative
abstract type DiffPeriodicBSplineBasis{T,K,D} <: CompactTranslationDict{T}
end

# For the B spline with degree 1 (hat functions) the MidpointEquispacedGrid does not lead to evaluation_matrix that is non singular
compatible_grid(b::DiffPeriodicBSplineBasis, grid::MidpointEquispacedGrid) = iseven(degree(b)) &&
    suppport(b)≈support(grid) && (length(b)==length(grid))

compatible_grid(b::DiffPeriodicBSplineBasis, grid::PeriodicEquispacedGrid) = isodd(degree(b)) &&
    suppport(b)≈support(grid) && (length(b)==length(grid))
    # we use a PeriodicEquispacedGrid in stead

interpolation_grid(b::DiffPeriodicBSplineBasis) = isodd(degree(b)) ?
    PeriodicEquispacedGrid(length(b), support(b)) :
    MidpointEquispacedGrid(length(b), support(b))

kernel_span(b::DiffPeriodicBSplineBasis) = Interval(zero(domaintype(b)), step(b)*convert(domaintype(b), degree(b)+1))

eval_kernel(b::DICT where DICT<:DiffPeriodicBSplineBasis{T,K,D}, x) where {T,K,D} =
    (n = length(b); sqrt(T(n))*n^D*evaluate_periodic_Bspline_derivative(Val{K}(), Val{D}(), n*x, n, T))

degree(b::DICT where DICT <: DiffPeriodicBSplineBasis{T,K} ) where {T,K} = K

Bdiff(b::DICT where DICT <: DiffPeriodicBSplineBasis{T,K,D} ) where {T,K,D} = D


abstract type PeriodicBSplineBasis{T,K} <: DiffPeriodicBSplineBasis{T,K,0}
end



function extension_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
    @assert degree(s1) == degree(s2)
    bspline_extension_operator(s1, s2; options...)
end

function restriction_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
    @assert degree(s1) == degree(s2)
    bspline_restriction_operator(s1, s2; options...)
end

function bspline_extension_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
    @assert 2*length(s1) == length(s2)
    _binomial_circulant(s2)*IndexExtensionOperator(s1, s2, 1:2:length(s2))
end

# The calculation done in this function is equivalent to finding the pseudoinverse of the bspline_extension_operator.
function bspline_restriction_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
    @assert length(s1) == 2*length(s2)
    r = _binomial_circulant(s1)
    e = eigvals(r)
    n = length(e)
    d = similar(e)
    eabs = map(abs, e)
    for i in 1:n>>1
        a = 2*(eabs[i]^2)/(eabs[i+n>>1]^2+eabs[i]^2)
        d[i] = a
        d[i+n>>1] = (2-a)
    end
    d .= d ./ e
    d[map(isnan,d)] .= 0
    IndexRestrictionOperator(s1,s2,1:2:length(s1))*CirculantOperator(real.(ifft(d)), s1, s1)
end

"""
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct BSplineTranslatesBasis{T,K,SCALED} <: PeriodicBSplineBasis{T,K}
    n               :: Int
end


Base.length(b::BSplineTranslatesBasis) = b.n
Base.size(b::BSplineTranslatesBasis) = (length(b),)

Base.similar(d::BSplineTranslatesBasis{S,K,SCALED}, ::Type{T}, n::Int) where {S,T,K,SCALED} = BSplineTranslatesBasis{T,K,SCALED}(n)

BSplineTranslatesBasis(n::Int, degree::Int, ::Type{T} = Float64; options...) where {T} =
    BSplineTranslatesBasis{T}(n, degree; options...)

BSplineTranslatesBasis{T}(n::Int, degree::Int; options...) where {T} = BSplineTranslatesBasis{T,degree}(n; options...)
BSplineTranslatesBasis{T,degree}(n::Int; scaled=true) where {T,degree} = BSplineTranslatesBasis{T,degree,scaled}(n)



eval_kernel(b::BSplineTranslatesBasis{T,K,true}, x) where {K,T} =
    (n = length(b); sqrt(T(n))*evaluate_periodic_Bspline(Val{K}(), n*x, n, T))
eval_kernel(b::BSplineTranslatesBasis{T,K,false}, x) where {K,T} =
    (n = length(b); evaluate_periodic_Bspline(Val{K}(), n*x, n, T))

scaled(b::BSplineTranslatesBasis{T,K,SCALED}) where {T,K,SCALED} = SCALED

name(b::BSplineTranslatesBasis) = name(typeof(b))*" (B spline of degree $(degree(b)))"

==(b1::BSplineTranslatesBasis{T1,K1}, b2::BSplineTranslatesBasis{T2,K2}) where {K1,K2,T1,T2} = T1==T2 && K1==K2 && length(b1)==length(b2)

instantiate(::Type{BSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = BSplineTranslatesBasis{T}(n,3)

resize(b::BSplineTranslatesBasis{T}, n::Int) where {T} = BSplineTranslatesBasis{T,degree(b)}(n)

innerproduct_native(d1::BSplineTranslatesBasis, i, d2::BSplineTranslatesBasis, j, measure::FourierMeasure; warnslow=false, options...)  =
     BasisFunctions.default_dict_innerproduct(d1, i, d2, j, measure; warnslow=warnslow, options...)

function _binomial_circulant(s::BSplineTranslatesBasis{T,K,SCALED}) where {T,K,SCALED}
    A = coefficienttype(s)
    c = zeros(A, length(s))
    for k in 1:K+2
        c[k] = binomial(K+1, k-1)
    end
    if SCALED
        sqrt(A(1//2))/(1<<K)*CirculantOperator(c, s)
    else
        A(1)/(1<<K)*CirculantOperator(c, s)
    end
end

function firstgramcolumn(dict::BSplineTranslatesBasis, measure::FourierMeasure; T=coefficienttype(s, domaintype(measure)), options...)
    firstcolumn = zeros(T, length(dict))
    for i in 1:length(dict)
        firstcolumn[i] = firstgramcolumnelement(dict, measure, i; T=T, options...)
    end
    firstcolumn
end

iscompatible(dict::BSplineTranslatesBasis, grid::MidpointEquispacedGrid) = iseven(degree(dict)) && isperiodic_compatible_grid(dict, grid)
iscompatible(dict::BSplineTranslatesBasis, grid::PeriodicEquispacedGrid) = isodd(degree(dict)) && isperiodic_compatible_grid(dict, grid)

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

Base.length(b::DiffBSplineTranslatesBasis) = b.n
Base.size(b::DiffBSplineTranslatesBasis) = (length(b),)

Base.similar(d::DiffBSplineTranslatesBasis{S,K,D}, ::Type{T}, n::Int) where {S,T,K,D} =
    DiffBSplineTranslatesBasis{T,K,D}(n)



DiffBSplineTranslatesBasis(n::Int, degree::Int, diff::Int, ::Type{T} = Float64; options...) where {T} =
    DiffBSplineTranslatesBasis{T}(n, degree, diff; options...)

DiffBSplineTranslatesBasis{T}(n::Int, degree::Int, diff::Int; options...) where {T} = DiffBSplineTranslatesBasis{T,degree,diff}(n; options...)
DiffBSplineTranslatesBasis{T,degree}(n::Int; D=0) where {T,degree} = DiffBSplineTranslatesBasis{T,degree,D}(n)

instantiate(::Type{DiffBSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = DiffBSplineTranslatesBasis(n,3,1,T)

resize(b::DiffBSplineTranslatesBasis, n::Int) = DiffBSplineTranslatesBasis(n, degree(b), Bdiff(b), codomaintype(b))
