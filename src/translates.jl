

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
function kernel_span(::CompactTranslationDict) end

length(dict::CompactTranslationDict) = dict.n

is_biorthogonal(::CompactTranslationDict) = true

is_basis(::CompactTranslationDict) = true

name(b::CompactTranslationDict) = "Dictionary of equispaced translates of a kernel function"
name(::Type{B}) where {B<:CompactTranslationDict}= "Dictionary of equispaced translates of a kernel function"

# Indices of translates naturally range from 0 to n-1
const TransIndex = ShiftedIndex{1}
ordering(b::CompactTranslationDict) = ShiftedIndexList{1}(length(b))

has_unitary_transform(::CompactTranslationDict) = false

support(::CompactTranslationDict{T}) where T = UnitInterval{T}()

support(b::CompactTranslationDict, idx) = support(b, native_index(b, idx))

support(b::CompactTranslationDict{T}, idx::TransIndex) where {T} = periodize_interval(value(idx)*stepsize(b)+kernel_span(b), support(b))

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
        (interval(infimum(outdomain), supremum(indomain)-per))∪(interval(infimum(indomain), supremum(outdomain)))
    end
end

period(::CompactTranslationDict{T}) where T = T(1)

stepsize(dict::CompactTranslationDict) = period(dict)/length(dict)

unsafe_eval_element(b::CompactTranslationDict, idxn::TransIndex, x::Real) =
    eval_kernel(b, x-idxn*stepsize(b))

function overlapping_elements(b::CompactTranslationDict, x::Real)
   indices = ceil(Int, (x-supremum(support(b)))/stepsize(b)):floor(Int, (x-infimum(support(b)))/stepsize(b))
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

has_grid(::CompactTranslationDict) = true

grid(dict::CompactTranslationDict) = PeriodicEquispacedGrid(length(dict), support(dict))

has_grid_transform(b::CompactTranslationDict, gb, grid::AbstractEquispacedGrid) =
    compatible_grid(b, grid)

compatible_grid(b::CompactTranslationDict, grid::AbstractEquispacedGrid) =
    periodic_compatible_grid(b, grid)

isdyadic(n::Integer) = (n == 2^(ndyadicscales(n)))
ndyadicscales(n::Integer) = round(Int, log2(n))
function periodic_compatible_grid(b::Dictionary, grid::AbstractEquispacedGrid)
    l1 = length(b)
    l2 = length(grid)
    l1 > l2 && ((l2,l1) = (l1, l2))
    n = l2/l1
    nInt = round(Int, n)
    (1+(infimum(support(b)) - leftendpoint(grid))≈1) && (1+(supremum(support(b)) - rightendpoint(grid))≈1) && isdyadic(nInt) && (n≈nInt)
end

approx_length(b::CompactTranslationDict, n::Int) = ceil(Int,n/length(b))*length(b)

native_nodes(b::CompactTranslationDict) = [k*stepsize(b) for k in 0:length(b)]

transform_from_grid(src, dest::CompactTranslationDict, grid; options...) =
    inv(transform_to_grid(dest, src, grid; options...))

function transform_to_grid(src::CompactTranslationDict, dest, grid; options...)
    @assert compatible_grid(src, grid)
    CirculantOperator(src, dest, sample(grid, x->eval_kernel(s, x)); options...)
end

function grid_evaluation_operator(s::CompactTranslationDict, dgs::GridBasis, grid::AbstractEquispacedGrid; TYPE=IndexableVerticalBandedOperator, options...)
    # r = nothing
    # if TYPE != IndexableVerticalBandedOperator
    #     if periodic_compatible_grid(s, grid)
    #         lg = length(grid)
    #         ls = length(s)
    #         sampling_factor, rem = divrem(lg, ls)
    #         @assert rem == 0
    #         if lg == ls
    #             r = CirculantOperator(s, dgs, sample(grid, x->eval_kernel(s, x)); options...)
    #         elseif lg > ls
    #             r = CirculantOperator(dgs, dgs, sample(grid, x->eval_kernel(s, x)); options...)*IndexExtensionOperator(s, dgs, 1:sampling_factor:length(dgs))
    #         elseif lg < ls && has_extension(grid)
    #             r = IndexRestrictionOperator(s, dgs, 1:sampling_factor:length(s))*CirculantOperator(s, s, sample(extend(grid, sampling_factor), x->eval_kernel(s, x)); options...)
    #         else
    #             r = default_evaluation_operator(s, dgs; options...)
    #         end
    #     else
    #         warn("slow evaluation operator")
    #         r = default_evaluation_operator(s, dgs; options...)
    #     end
    #     if TYPE == SparseOperator
    #         return SparseOperator(r; options...)
    #     else
    #         return r
    #     end
    # else
    lg = length(grid)
    ls = length(s)
    sampling_factor, rem = divrem(lg, ls)
    if rem == 0
        firstcolumn = sample(grid, x->eval_kernel(s, x))
        a, offset = _get_array_offset(firstcolumn)
        IndexableVerticalBandedOperator(s, dgs, a, sampling_factor, offset-1)
    else
        warn("slow evaluation operator")
        default_evaluation_operator(s, dgs; options...)
    end
end

function _get_array_offset(a)
    b = a.!=0
    f = findfirst(b)
    if f==1
        if b[end]
            f = findlast(.!b)+1
            L = sum(b)
            vcat(a[f:end],a[1:L-length(a)+f]), f
        else
            a[f:f+sum(b)-1], f
        end
    else
        a[f:f+sum(b)-1], f
    end
end

function grid_evaluation_operator(s::D, dgs::GridBasis, grid::ProductGrid;
        options...) where {D<: TensorProductDict{N,DT,S,T} where {N,DT <: NTuple{N,CompactTranslationDict} where N,S,T}}
    tensorproduct([grid_evaluation_operator(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)
end

Gram(s::CompactTranslationDict; options...) = CirculantOperator(s, s, primalgramcolumn(s; options...))

function UnNormalizedGram(s::CompactTranslationDict, oversampling = 1)
    grid = BasisFunctions.oversampled_grid(s, oversampling)
    CirculantOperator(evaluation_operator(s, grid)'*evaluation_operator(s, grid))
end

# All inner products between elements of CompactTranslationDict are known by the first column of the (circulant) gram matrix.
function primalgramcolumn(s::CompactTranslationDict; options...)
    n = length(s)
    result = zeros(coeftype(s), n)
    for i in 1:length(result)
        result[i] = primalgramcolumnelement(s, i; options...)
    end
    result
end

primalgramcolumnelement(s::CompactTranslationDict, i::Int; options...) =
    defaultprimalgramcolumnelement(s, i; options...)

defaultprimalgramcolumnelement(s::Dictionary1d, i::Int; options...)  = dot(s, 1, i; options...)

"""
  The function set that correspands to the dual set ``Ψ={ψ_i}_{i∈ℕ}`` of the given set ``Φ={ϕ_i}_{i∈ℕ}``.

"""
dual(set::CompactTranslationDict; options...) =
    DualCompactTranslationDict(set; options...)
