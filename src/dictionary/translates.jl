
"""
Dictionary periodic on [0,1] consisting of equispaced translates of one generating function.
"""
abstract type CompactTranslationDict{T} <: Dictionary{T,T}
end

"""
Evaluate the generating function of in x.
"""
function eval_kernel(::CompactTranslationDict, x) end

"""
The span of the kernel.
"""
function kernel_support(::CompactTranslationDict) end

isbasis(::CompactTranslationDict) = true

name(b::CompactTranslationDict) = "Dictionary of equispaced translates of a kernel function"
name(::Type{B}) where {B<:CompactTranslationDict}= "Dictionary of equispaced translates of a kernel function"

# Indices of translates naturally range from 0 to n-1
const TransIndex = ShiftedIndex{1}
ordering(b::CompactTranslationDict) = ShiftedIndexList{1}(length(b))

support(::CompactTranslationDict{T}) where T = UnitInterval{T}()

support(b::CompactTranslationDict, idx) = support(b, native_index(b, idx))

support(b::CompactTranslationDict{T}, idx::TransIndex) where {T} = periodize_interval(value(idx)*step(b)+kernel_support(b), support(b))

function periodize_interval(indomain::AbstractInterval, outdomain::AbstractInterval)
    per = (supremum(outdomain) - infimum(outdomain))
    if ((supremum(indomain) - infimum(indomain)) >= per)
        return outdomain
    end
    if (infimum(indomain) ∉ outdomain)
        error("Case not yet implemented.")
    end
    if supremum(indomain) ∈ outdomain
        indomain
    else
        UnionDomain((Interval(infimum(outdomain), supremum(indomain)-per)),(Interval(infimum(indomain), supremum(outdomain))))
    end
end

period(::CompactTranslationDict{T}) where T = T(1)

step(dict::CompactTranslationDict) = period(dict)/length(dict)

unsafe_eval_element(b::CompactTranslationDict, idxn::TransIndex, x::Real) =
    eval_kernel(b, x-idxn*step(b))

function overlapping_elements(b::CompactTranslationDict, x::Real)
   indices = ceil(Int, (x-supremum(support(b)))/step(b)):floor(Int, (x-infimum(support(b)))/step(b))
   Set(mod(i, length(b))+1 for i in indices)
end

function eval_expansion(b::CompactTranslationDict, coef, x)
    z = zero(typeof(x))
    for idx = overlapping_elements(b, x)
        idxn = native_index(b, idx)
        z = z + coef[idx] * unsafe_eval_element1(b, idxn, x)
    end
    z
end

hasinterpolationgrid(::CompactTranslationDict) = true

interpolation_grid(dict::CompactTranslationDict) = PeriodicEquispacedGrid(length(dict), support(dict))
function oversampling_grid(dict::CompactTranslationDict, L)
    L = ceil(Int, L/length(dict))*length(dict)
    interpolation_grid(resize(dict, L))
end

hasgrid_transform(b::CompactTranslationDict, gb, grid::AbstractEquispacedGrid) =
    compatible_translationgrid(b, grid)

iscompatible(b::CompactTranslationDict, grid::AbstractEquispacedGrid) =
    isperiodic_compatible_translationgrid(b, grid)

isdyadic(n::Integer) = (n == 2^(ndyadicscales(n)))
ndyadicscales(n::Integer) = round(Int, log2(n))
isperiodic_compatible_translationgrid(b::Dictionary, grid::AbstractGrid) = false
function isperiodic_compatible_translationgrid(b::CompactTranslationDict, grid::AbstractEquispacedGrid)
    l1 = length(b)
    l2 = length(grid)
    l1 > l2 && ((l2,l1) = (l1, l2))
    n = l2/l1
    nInt = round(Int, n)
    support(b)≈support(grid) && (n≈nInt)
end

approx_length(b::CompactTranslationDict, n::Int) = ceil(Int,n/length(b))*length(b)

transform_from_grid(src, dest::CompactTranslationDict, grid; options...) =
    inv(transform_to_grid(dest, src, grid; options...))

function transform_to_grid(src::CompactTranslationDict, dest, grid; options...)
    @assert hasgrid_transform(src, dest, grid)
    CirculantOperator(src, dest, sample(grid, x->eval_kernel(src, x)); options...)
end

function grid_evaluation_operator(s::CompactTranslationDict, dgs::GridBasis, grid::AbstractEquispacedGrid; T=op_eltype(s, dgs), TYPE=VerticalBandedOperator, warnslow = BasisFunctions.BF_WARNSLOW, options...)
    lg = length(grid)
    ls = length(s)
    sampling_factor, rem = divrem(lg, ls)
    if rem == 0
        firstcolumn = sample(grid, x->eval_kernel(s, x))
        a, offset = _get_array_offset(firstcolumn)
        BasisFunctions.VerticalBandedOperator(s, dgs, a, sampling_factor, offset-1; T=T)
    else
        warnslow && (@warn("slow evaluation operator"))
        dense_evaluation_operator(s, dgs; options...)
    end
end

function _get_array_offset(a)
    b = a.!=0
    f = findfirst(b)
    (f==nothing) && (f=0)

    if f==1
        if b[end]
            f = findlast(.!b)
            (f == nothing) ? (return (a, 1)) : f += 1
            L = sum(b)
            vcat(a[f:end],a[1:L-length(a)+f]), f
        else
            a[f:f+sum(b)-1], f
        end
    else
        a[f:f+sum(b)-1], f
    end
end

hasmeasure(dict::CompactTranslationDict) = true
measure(b::CompactTranslationDict{T}) where T = FourierMeasure{T}()

gramoperator(dict::CompactTranslationDict, measure::Union{GenericLebesgueMeasure,LegendreMeasure,FourierMeasure}; options...) =
    _translatescirculantoperator(dict, measure)

function gramoperator(dict::CompactTranslationDict, measure::BasisFunctions.DiscreteMeasure, grid::AbstractEquispacedGrid, weights::FillArrays.AbstractFill; options...)
    if isperiodic_compatible_translationgrid(dict, grid)
        CirculantOperator(BasisFunctions.default_mixedgramoperator_discretemeasure(dict, dict, measure, grid, weights; options...))
    else
        BasisFunctions.default_mixedgramoperator_discretemeasure(dict, dict, measure, grid, weights; options...)
    end
end

_translatescirculantoperator(dict::Dictionary, measure::Measure; T = coefficienttype(dict), options...) =
	CirculantOperator(firstgramcolumn(dict, measure; T=T, options...), dict, dict; T=T)

function firstgramcolumn(dict::Dictionary, measure::Measure; T = coefficienttype(dict), options...)
    firstcolumn = zeros(T, length(dict))
    for (index,i) in enumerate(ordering(dict))
        firstcolumn[index] = innerproduct(dict, i, dict, ordering(dict)[1], measure; options...)
    end
    firstcolumn
end
