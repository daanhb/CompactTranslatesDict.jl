const CompactTranslationDict = PeriodicEquispacedTranslates

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

hasgrid_transform(dict::PeriodicEquispacedTranslates, _, grid::AbstractEquispacedGrid) =
    size(dict)==size(grid) && compatible_translationgrid(typeof(dict), grid)

transform_from_grid(src, dest::PeriodicEquispacedTranslates, grid; options...) =
    inv(transform_to_grid(dest, src, grid; options...))

function transform_to_grid(src::PeriodicEquispacedTranslates, dest, grid::AbstractEquispacedGrid; options...)
    @assert hasgrid_transform(src, dest, grid)
    CirculantOperator(src, dest, sample(grid, x->eval_kernel(src, x)); options...)
end
using BasisFunctions: VerticalBandedMatrix

function grid_evaluation_operator(s::PeriodicEquispacedTranslates, dgs::GridBasis, grid::AbstractEquispacedGrid;
            T=op_eltype(s, dgs), options...)
    lg = length(grid)
    ls = length(s)
    sampling_factor, rem = divrem(lg, ls)
    if rem == 0
        firstcolumn = sample(grid, x->(unsafe_eval_element(s, 1, x)))
        firstcolumn = convert.(T, firstcolumn)
        a, offset = _get_array_offset(firstcolumn)
        M = VerticalBandedMatrix{T}(length(dgs), length(s), a, sampling_factor, offset-1)
        ArrayOperator(M, s, dgs)
        # VerticalBandedOperator(s, dgs, a, sampling_factor, offset-1; T=T)
    else
        @debug "slow evaluation operator"
        # Not type stable if we allow this branch
        BasisFunctions.dense_evaluation_operator(s, dgs; options...)
        # error("Not type stable if we allow this branch")
    end
end

gramoperator(dict::PeriodicEquispacedTranslates, measure::GenericLebesgueMeasure; options...) =
    _translatescirculantoperator(dict, measure)

if isdefined(BasisFunctions, :isgramcompatible)
    import BasisFunctions: isgramcompatible
else
    isgramcompatible(dict::Dictionary, grid::AbstractGrid) = false
end

function isgramcompatible(b::PeriodicEquispacedTranslates, grid::AbstractEquispacedGrid)
    l1 = length(b)
    l2 = length(grid)
    l1 > l2 && ((l2,l1) = (l1, l2))
    n = l2/l1
    nInt = round(Int, n)
    support(b)≈support(grid) && (n≈nInt)
end

function gramoperator(dict::PeriodicEquispacedTranslates, measure::DiscreteMeasure, grid::AbstractEquispacedGrid, weights::FillArrays.AbstractFill;
        options...)
    if isgramcompatible(dict, grid)
        CirculantOperator(default_mixedgramoperator_discretemeasure(dict, dict, measure, grid, weights; options...))
    else
        default_mixedgramoperator_discretemeasure(dict, dict, measure, grid, weights; options...)
    end
end


function gramoperator(dict::PeriodicEquispacedTranslates, measure::Union{GenericLebesgueMeasure,LegendreMeasure,FourierMeasure};
        T = coefficienttype(dict), options...)
    @assert support(dict) ≈ support(measure)
    CirculantOperator(firstgramcolumn(dict, measure; T=T, options...), dict, dict; T=T)
end

function firstgramcolumn(dict::Dictionary, measure::Measure; T = coefficienttype(dict), options...)
    firstcolumn = zeros(T, length(dict))
    for (index,i) in enumerate(ordering(dict))
        firstcolumn[index] = innerproduct(dict, i, dict, ordering(dict)[1], measure; options...)
    end
    firstcolumn
end

# The evaluation of the basis functions explicitly periodizes the kernel function
# by summing over its translates.
function unsafe_eval_element(dict::PeriodicEquispacedTranslates{T,S,:sum}, idx, x) where {T,S}
    c = translationgrid(dict)[idx]
    per = period(dict)
    A,B = extrema(kernel_support(dict))

    z = eval_kernel(dict, x - c)
	# Now evaluate the periodic extension. We add and subtract the period
	# repeatedly until we are beyond the support of the kernel.
    x1 = x + per
    while (x1 <= c + B)
        z += eval_kernel(dict, x1-c)
        x1 += per
    end
    x2 = x - per
    while (x2 >= c + A)
        z += eval_kernel(dict, x2-c)
		x2 -= per
    end
    z
end

unsafe_eval_element(dict::PeriodicEquispacedTranslates{T,S,PERIODIZATION}) where {T,S,PERIODIZATION} =
    error("Periodization $PERIODIZATION not known.")

function LinearAlgebra.norm(dict::PeriodicEquispacedTranslates{T,S,:norm}, x, y) where {T,S}
    per = period(dict)
    ELT = typeof(per)
    # r = max(abs.(extrema(kernel_support(dict)))...)
    # t = r^2/(1-cos(2ELT(pi)/per*r))
    sqrt((1-cos(2ELT(pi)/per*(x - y))))
end
function unsafe_eval_element(dict::PeriodicEquispacedTranslates{T,S,:norm}, idx, x) where {T,S}
    c = translationgrid(dict)[idx]
    z = eval_kernel(dict, norm(dict, x, c))
end
