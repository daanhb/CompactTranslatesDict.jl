# translates_of_bsplines.jl
using CardinalBSplines
abstract type DiffPeriodicBSplineBasis{K,T,D} <: CompactTranslationDict{T}
end

# For the B spline with degree 1 (hat functions) the MidpointEquispacedGrid does not lead to evaluation_matrix that is non singular
compatible_grid(b::DiffPeriodicBSplineBasis, grid::MidpointEquispacedGrid) = iseven(degree(b)) &&
    (1+(infimum(support(b)) - leftendpoint(grid))≈1) && (1+(supremum(support(b)) - rightendpoint(grid))≈1) && (length(b)==length(grid))

compatible_grid(b::DiffPeriodicBSplineBasis, grid::PeriodicEquispacedGrid) = isodd(degree(b)) &&
    (1+(infimum(support(b)) - leftendpoint(grid))≈1) && (1+(supremum(support(b)) - rightendpoint(grid))≈1) && (length(b)==length(grid))
    # we use a PeriodicEquispacedGrid in stead

grid(b::DiffPeriodicBSplineBasis) = isodd(degree(b)) ? PeriodicEquispacedGrid(length(b), support(b)) : MidpointEquispacedGrid(length(b), support(b))

kernel_span(b::DiffPeriodicBSplineBasis) = Interval(domaintype(b)(0), stepsize(b)*domaintype(b)(degree(b)+1))

eval_kernel(b::DICT where DICT<:DiffPeriodicBSplineBasis{K,T,D}, x) where {K,T,D} = (n = length(b); sqrt(T(n))*n^D*diff_evaluate_periodic_Bspline(Val{K}, Val{D}, n*x, n, T))::T


abstract type PeriodicBSplineBasis{K,T} <: DiffPeriodicBSplineBasis{K,T,0}
end

degree(b::DICT where DICT <: DiffPeriodicBSplineBasis{K} ) where {K} = K

Bdiff(b::DICT where DICT <: DiffPeriodicBSplineBasis{K,T,D} ) where {K,T,D} = D

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
    e = eigenvalues(r)
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

    IndexRestrictionOperator(s1,s2,1:2:length(s1))*CirculantOperator(s1, s1, DiagonalOperator(d))
end

"""
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct BSplineTranslatesBasis{K,T,SCALED} <: PeriodicBSplineBasis{K,T}
    n               :: Int
end

Base.length(b::BSplineTranslatesBasis) = b.n
Base.size(b::BSplineTranslatesBasis) = (length(b),)

Base.similar(d::BSplineTranslatesBasis, ::Type{T}, n::Int) where T = BSplineTranslatesBasis(n, degree(d), T)

BSplineTranslatesBasis(n::Int, DEGREE::Int, ::Type{T} = Float64; scaled = true) where {T} = (scaled) ?
    BSplineTranslatesBasis{DEGREE,T,true}(n) :
    BSplineTranslatesBasis{DEGREE,T,false}(n)

eval_kernel(b::BSplineTranslatesBasis{K,T,true}, x) where {K,T} = (n = length(b); sqrt(T(n))*evaluate_periodic_Bspline(Val{K}, n*x, n, T))::T
eval_kernel(b::BSplineTranslatesBasis{K,T,false}, x) where {K,T} = (n = length(b); evaluate_periodic_Bspline(Val{K}, n*x, n, T))::T

scaled(b::BSplineTranslatesBasis{K,T,SCALED}) where {K,T,SCALED} = SCALED

name(b::BSplineTranslatesBasis) = name(typeof(b))*" (B spline of degree $(degree(b)))"

==(b1::BSplineTranslatesBasis{K1,T1}, b2::BSplineTranslatesBasis{K2,T2}) where {K1,K2,T1,T2} = T1==T2 && K1==K2 && length(b1)==length(b2)

instantiate(::Type{BSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = BSplineTranslatesBasis(n,3,T)

resize(b::BSplineTranslatesBasis{K,T}, n::Int) where {K,T} = BSplineTranslatesBasis(n, degree(b), T)

function _binomial_circulant(s::BSplineTranslatesBasis{K,T,SCALED}) where {K,T,SCALED}
    A = coefficienttype(s)
    c = zeros(A, length(s))
    for k in 1:K+2
        c[k] = binomial(K+1, k-1)
    end
    if SCALED
        sqrt(A(1//2))/(1<<K)*CirculantOperator(s, c)
    else
        A(1)/(1<<K)*CirculantOperator(s, c)
    end
end

function primalgramcolumnelement(s::BSplineTranslatesBasis{K,T,SCALED}, i::Int; options...) where {K,T,SCALED}
    r = 0
    A = coefficienttype(s)
    # If size of functionspace is too small there is overlap and we can not use the
    # function squared_spline_integral which assumes no overlap.
    # Use integration as long as there is no more efficient way is implemented.
    if length(s) <= 2K+1
        return defaultprimalgramcolumnelement(s, i; options...)
    else
        # squared_spline_integral gives the exact integral (in a rational number)
        if i==1
            r = CardinalBSplines.squared_spline_integral(K)
        elseif 1 < i <= K+1
            r = CardinalBSplines.shifted_spline_integral(K,i-1)
        elseif i > length(s)-K
            r = CardinalBSplines.shifted_spline_integral(K,length(s)-i+1)
        end
    end
    if SCALED
        A(r)
    else
        A(r)/length(s)
    end
end

"""
  Basis consisting of differentiated dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct DiffBSplineTranslatesBasis{K,T,D} <: DiffPeriodicBSplineBasis{K,T,D}
    n               :: Int
end

Base.length(b::DiffBSplineTranslatesBasis) = b.n
Base.size(b::DiffBSplineTranslatesBasis) = (length(b),)

Base.similar(d::DiffBSplineTranslatesBasis, ::Type{T}, n::Int) where T = DiffBSplineTranslatesBasis(N, degree(d), Bdiff(d), T)

DiffBSplineTranslatesBasis(n::Int, DEGREE::Int, D::Int, ::Type{T} = Float64) where {T} =
    DiffBSplineTranslatesBasis{DEGREE,T,D}(n)

instantiate(::Type{DiffBSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = DiffBSplineTranslatesBasis(n,3,1,T)

resize(b::DiffBSplineTranslatesBasis, n::Int) = DiffBSplineTranslatesBasis(n, degree(b), Bdiff(b), codomaintype(b))
