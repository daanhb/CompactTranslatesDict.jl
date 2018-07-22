# translates_of_bsplines.jl
using CardinalBSplines
abstract type PeriodicBSplineBasis{K,T} <: CompactTranslationDict{T}
end

degree(b::B) where {K,T, B<:PeriodicBSplineBasis{K,T}} = K

BasisFunctions.Gram(b::PeriodicBSplineBasis; options...) = CirculantOperator(b, b, primalgramcolumn(b; options...); options...)

function BasisFunctions.extension_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
    @assert degree(s1) == degree(s2)
    bspline_extension_operator(s1, s2; options...)
end

function BasisFunctions.restriction_operator(s1::PeriodicBSplineBasis, s2::PeriodicBSplineBasis; options...)
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
    e = BasisFunctions.eigenvalues(r)
    n = length(e)
    d = similar(e)
    eabs = map(abs, e)
    for i in 1:n>>1
      a = 2*(eabs[i]^2)/(eabs[i+n>>1]^2+eabs[i]^2)
      d[i] = a
      d[i+n>>1] = (2-a)
    end
    d = d ./ e
    d[map(isnan,d)] .= 0

    IndexRestrictionOperator(s1,s2,1:2:length(s1))*CirculantOperator(s1, s1, DiagonalOperator(d))
end

"""
  Basis consisting of dilated, translated, and periodized cardinal B splines on the interval [0,1].
"""
struct BSplineTranslatesBasis{K,T,SCALED} <: PeriodicBSplineBasis{K,T}
  n               :: Int
end

BSplineTranslatesBasis(n::Int, DEGREE::Int, ::Type{T} = Float64; scaled = false) where {T} = (scaled) ?
    BSplineTranslatesBasis{DEGREE,T,true}(n) :
    BSplineTranslatesBasis{DEGREE,T,false}(n)

eval_kernel(b::BSplineTranslatesBasis{K,T,true}, x) where {K,T} = (n = length(b); sqrt(n)*evaluate_periodic_Bspline(Val{K}, n*x, n, T))
eval_kernel(b::BSplineTranslatesBasis{K,T,false}, x) where {K,T} = (n = length(b); evaluate_periodic_Bspline(Val{K}, n*x, n, T))

kernel_span(b::BSplineTranslatesBasis{K,T}) where {K,T} = interval(domaintype(b)(0), stepsize(b)*real(T)(degree(b)+1))

name(b::BSplineTranslatesBasis) = name(typeof(b))*" (B spline of degree $(degree(b)))"

==(b1::BSplineTranslatesBasis{K1,T1}, b2::BSplineTranslatesBasis{K2,T2}) where {K1,K2,T1,T2} = T1==T2 && K1==K2 && length(b1)==length(b2)

BasisFunctions.instantiate(::Type{BSplineTranslatesBasis}, n::Int, ::Type{T}) where {T} = BSplineTranslatesBasis(n,3,T)

BasisFunctions.dict_promote_domaintype(b::BSplineTranslatesBasis{K,T}, ::Type{S}) where {K,S,T} = BSplineTranslatesBasis(length(b),K, S)

BasisFunctions.resize(b::BSplineTranslatesBasis{K,T}, n::Int) where {K,T} = BSplineTranslatesBasis(n, degree(b), T)

# For the B spline with degree 1 (hat functions) the MidpointEquispacedGrid does not lead to evaluation_matrix that is non singular
BasisFunctions.compatible_grid(b::BSplineTranslatesBasis{K}, grid::MidpointEquispacedGrid) where {K} = iseven(K) &&
    (1+(infimum(support(b)) - leftendpoint(grid))≈1) && (1+(supremum(support(b)) - rightendpoint(grid))≈1) && (length(b)==length(grid))
BasisFunctions.compatible_grid(b::BSplineTranslatesBasis{K}, grid::PeriodicEquispacedGrid) where {K} = isodd(K) &&
    (1+(infimum(support(b)) - leftendpoint(grid))≈1) && (1+(supremum(support(b)) - rightendpoint(grid))≈1) && (length(b)==length(grid))
    # we use a PeriodicEquispacedGrid in stead
BasisFunctions.grid(b::BSplineTranslatesBasis{K}) where {K} = isodd(K) ? PeriodicEquispacedGrid(length(b), support(b)) : MidpointEquispacedGrid(length(b), support(b))

function _binomial_circulant(s::BSplineTranslatesBasis{K,T,SCALED}) where {K,T,SCALED}
    A = coeftype(s)
    c = zeros(A, length(s))
    for k in 1:K+2
        c[k] = binomial(K+1, k-1)
    end
    if SCALED
        sqrt(A(2))/(1<<K)*CirculantOperator(s, c)
    else
        A(1)/(1<<K)*CirculantOperator(s, c)
    end
end

function primalgramcolumnelement(s::BSplineTranslatesBasis{K,T,SCALED}, i::Int; options...) where {K,T,SCALED}
    r = 0
    A = coeftype(s)
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
