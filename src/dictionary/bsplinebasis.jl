using CardinalBSplines: AbstractBSpline
using BasisFunctions.GridArrays: AbstractGrid
import BasisFunctions: support

struct BSplineBasis{T<:Real,K,SPLINE<:AbstractBSpline{K,T},GRID<:AbstractGrid} <: Translates{T,T}
    spline ::   SPLINE
    grid   ::   GRID

end

eval_kernel(dict::BSplineBasis, x) =
    dict.spline(x)
kernel_support(dict::BSplineBasis) =
    support(dict.spline)

support(spline::BSpline{K,S}) where {K,S} = Interval(zero(S), convert(S,K+1))

support(spline::CenteredBSpline) = Interval(-convert(S,K+1)/2, convert(S,K+1)/2)
