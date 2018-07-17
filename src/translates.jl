

"""
Dictionary consisting of translates of one generating function.
"""
abstract type TranslationDict{T} <: Dictionary{T,T}
end


BasisFunctions.length(set::TranslationDict) = set.n

BasisFunctions.is_biorthogonal(::TranslationDict) = true

BasisFunctions.is_basis(::TranslationDict) = true

BasisFunctions.name(b::TranslationDict) = "Set of translates of a function"
BasisFunctions.name(::Type{B}) where {B<:TranslationDict}= "Set of translates of a function"

fun(b::TranslationDict) = b.fun

# Indices of translates naturally range from 0 to n-1
const TransIndex = BasisFunctions.ShiftedIndex{1}
BasisFunctions.ordering(b::TranslationDict) = BasisFunctions.ShiftedIndexList{1}(length(b))

BasisFunctions.has_unitary_transform(::TranslationDict) = false


"""
  Set consisting of n equispaced translates of a periodic function.

  The set can be written as ``\left\{T_k f\right\}_{k=0}^n``, where ``T_k f(x) = f(x-p/n)``.
  ``p`` is the period of the set, ``n`` is the number of elements.
"""
abstract type PeriodicTranslationDict{T} <: TranslationDict{T}
end



BasisFunctions.support(set::PeriodicTranslationDict{T}) where {T} = interval(set.a,set.b)

BasisFunctions.support(set::PeriodicTranslationDict, j::TransIndex) = support(set)


BasisFunctions.has_grid(::PeriodicTranslationDict) = true

BasisFunctions.grid(set::PeriodicTranslationDict) = PeriodicEquispacedGrid(length(set), support(set))

BasisFunctions.period(set::PeriodicTranslationDict) = supremum(support(set))-infimum(support(set))

BasisFunctions.stepsize(set::PeriodicTranslationDict) = BasisFunctions.period(set)/length(set)

BasisFunctions.has_grid_transform(b::PeriodicTranslationDict, gb, grid::AbstractEquispacedGrid) =
    compatible_grid(b, grid)

BasisFunctions.compatible_grid(b::PeriodicTranslationDict, grid::AbstractEquispacedGrid) =
    periodic_compatible_grid(b, grid)

BasisFunctions.approx_length(b::PeriodicTranslationDict, n::Int) = ceil(Int,n/length(b))*length(b)

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

BasisFunctions.native_nodes(b::PeriodicTranslationDict) = [k*stepsize(b) for k in 0:length(b)]

function BasisFunctions.transform_from_grid(src, dest::PeriodicTranslationDict, grid; options...)
    inv(transform_to_grid(dest, src, grid; options...))
end

function BasisFunctions.transform_to_grid(src::PeriodicTranslationDict, dest, grid; options...)
    @assert compatible_grid(src, grid)
    CirculantOperator(src, dest, sample(grid, fun(src)); options...)
end


function BasisFunctions.grid_evaluation_operator(s::PeriodicTranslationDict, dgs::GridBasis, grid::AbstractEquispacedGrid; sparse=true, options...)
    r = nothing
    if periodic_compatible_grid(s, grid)
        lg = length(grid)
        ls = length(s)
        if lg == ls
            r = CirculantOperator(s, dgs, sample(grid, fun(s)); options...)
        elseif lg > ls
            r = CirculantOperator(dgs, dgs, sample(grid, fun(s)); options...)*IndexExtensionOperator(s, dgs, 1:Int(lg/ls):length(dgs))
        elseif lg < ls && has_extension(grid)
            r = IndexRestrictionOperator(s, dgs, 1:Int(ls/lg):length(s))*CirculantOperator(s, s, sample(extend(grid, Int(ls/lg)), fun(s)); options...)
        else
            r = default_evaluation_operator(s, dgs; options...)
        end
    else
        r = default_evaluation_operator(s, dgs; options...)
    end
    if sparse
        return BasisFunctions.SparseOperator(r; options...)
    else
        return r
    end
end

function BasisFunctions.grid_evaluation_operator(s::D, dgs::GridBasis, grid::ProductGrid;
        options...) where {D<: TensorProductDict{N,DT,S,T} where {N,DT <: NTuple{N,PeriodicTranslationDict} where N,S,T}}
    tensorproduct([BasisFunctions.grid_evaluation_operator(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)
end

BasisFunctions.unsafe_eval_element(b::PeriodicTranslationDict, idxn::TransIndex, x::Real) =
    fun(b)(x-value(idxn)*stepsize(b))

eval_dualelement(b::PeriodicTranslationDict, idx::LinearIndex, x::Real) =
    eval_dualelement(b, native_index(b, idx), x)

eval_dualelement(b::PeriodicTranslationDict, idxn::TransIndex, x::Real) =
    eval_expansion(b, circshift(dualgramcolumn(b),value(idxn)), x)

BasisFunctions.Gram(s::PeriodicTranslationDict; options...) = CirculantOperator(s, s, primalgramcolumn(s; options...))

function BasisFunctions.UnNormalizedGram(s::PeriodicTranslationDict, oversampling = 1)
    grid = BasisFunctions.oversampled_grid(s, oversampling)
    BasisFunctions.CirculantOperator(evaluation_operator(s, grid)'*evaluation_operator(s, grid))
end

grammatrix(b::PeriodicTranslationDict; options...) = matrix(Gram(b; options...))

dualgrammatrix(b::PeriodicTranslationDict; options...) = matrix(inv(Gram(b; options...)))

# All inner products between elements of PeriodicTranslationDict are known by the first column of the (circulant) gram matrix.
function primalgramcolumn(s::PeriodicTranslationDict; options...)
    n = length(s)
    result = zeros(coeftype(s), n)
    for i in 1:length(result)
        result[i] = primalgramcolumnelement(s, i; options...)
    end
    result
end

primalgramcolumnelement(s::PeriodicTranslationDict, i::Int; options...) =
    defaultprimalgramcolumnelement(s, i; options...)

defaultprimalgramcolumnelement(s::Dictionary1d, i::Int; options...)  = dot(s, 1, i; options...)

function dualgramcolumn(s::PeriodicTranslationDict; options...)
    G = inv(Gram(s; options...))
    e = zeros(eltype(G),size(G,1))
    e[1] = 1
    real(G*e)
end

"""
  The function set that correspands to the dual set ``Ψ={ψ_i}_{i∈ℕ}`` of the given set ``Φ={ϕ_i}_{i∈ℕ}``.

"""
BasisFunctions.dual(set::PeriodicTranslationDict; options...) =
    DualPeriodicTranslationDict(set; options...)
"""
  The function set that correspands to the dual set of the given set.
"""
discrete_dual(set::PeriodicTranslationDict; options...) =
    DiscreteDualPeriodicTranslationDict(set; options...)

"""
  Set consisting of n translates of a compact and periodic function.

  The support of the function is [c_1,c_2], where c_1, c_2 ∈R, c_2-c_1 <= p, 0 ∈ [c_1,c_2],
  and p is the period of the function.
"""
abstract type CompactPeriodicTranslationDict{T} <: PeriodicTranslationDict{T}
end



"""
  Length of the function of a CompactPeriodicTranslationDict.
"""
length_compact_support{T}(b::CompactPeriodicTranslationDict{T})::real(T) = right_of_compact_function(b)-left_of_compact_function(b)

function left_of_compact_function end
function right_of_compact_function end

support_length_of_compact_function(f::CompactPeriodicTranslationDict) = right_of_compact_function(f::CompactPeriodicTranslationDict)-left_of_compact_function(f::CompactPeriodicTranslationDict)

function overlapping_elements(b::CompactPeriodicTranslationDict, x::Real)
   indices = ceil(Int, (x-CompactTranslatesDict.right_of_compact_function(b))/stepsize(b)):floor(Int, (x-CompactTranslatesDict.left_of_compact_function(b))/stepsize(b))
   Set(mod(i, length(b))+1 for i in indices)
end

BasisFunctions.support(b::CompactPeriodicTranslationDict, idx) = support(b, native_index(b, idx))

BasisFunctions.support(b::CompactPeriodicTranslationDict, idx::TransIndex) = value(idx)*BasisFunctions.stepsize(b)+compact_support(b)

compact_support(b::CompactPeriodicTranslationDict) = interval(left_of_compact_function(b),right_of_compact_function(b))

BasisFunctions.in_support(set::CompactPeriodicTranslationDict, idx::LinearIndex, x::Real) =
    in_support(set, native_index(set, idx), x)

BasisFunctions.in_support(set::CompactPeriodicTranslationDict, idx::TransIndex, x::Real) =
    in_compact_support(set, idx, x)

BasisFunctions.unsafe_eval_element(b::CompactPeriodicTranslationDict, idx::TransIndex, x::Real) =
    eval_compact_element(b, idx, x)

BasisFunctions.eval_expansion(b::CompactPeriodicTranslationDict, coef, x::Real) =
    eval_compact_expansion(b, coef, x)

function eval_compact_element(b::CompactPeriodicTranslationDict{T}, idx::TransIndex, x::Real) where {T}
    !in_support(b, idx, x) ? zero(T) : fun(b)(x-value(idx)*stepsize(b))
end

function in_compact_support(set::CompactPeriodicTranslationDict, idx::TransIndex, x::Real)
	per = BasisFunctions.period(set)
        A = in(x,support(set))
	B = in(x,support(set,idx)) || in(x-per,support(set,idx)) || in(x+per,support(set,idx))
	A && B
end

function eval_compact_expansion(b::CompactPeriodicTranslationDict, coef, x)
    z = zero(typeof(x))
    for idx = overlapping_elements(b, x)
        idxn = native_index(b, idx)
        z = z + coef[idx] * BasisFunctions.unsafe_eval_element(b, idxn, x)
    end
    z
end

"""
  Set of translates of a function f that is a linear combination of basis function of an other set of translates.

  `f(x) = ∑ coefficients(set)[k] * superdictionary(set)[k](x)`
"""
abstract type LinearCombinationOfPeriodicTranslationDict{PSoT<:PeriodicTranslationDict, T} <: PeriodicTranslationDict{T}
end



BasisFunctions.coefficients(b::LinearCombinationOfPeriodicTranslationDict) = b.coefficients

for op in (:length, :has_grid, :grid, :support)
    @eval BasisFunctions.$op(b::LinearCombinationOfPeriodicTranslationDict) = BasisFunctions.$op(superdict(b))
end

function fun(b::LinearCombinationOfPeriodicTranslationDict)
    x->eval_expansion(superdict(b), real(BasisFunctions.coefficients(b)), CardinalBSplines.periodize(x, BasisFunctions.period(superdict(b))))
end

==(b1::LinearCombinationOfPeriodicTranslationDict, b2::LinearCombinationOfPeriodicTranslationDict) =
    superdict(b1)==superdict(b2) && BasisFunctions.coefficients(b1) ≈ BasisFunctions.coefficients(b2)

change_of_basis(b::LinearCombinationOfPeriodicTranslationDict; options...) =
    wrap_operator(superdict(b), b, inv(change_of_basis(superdict(b), typeof(b))))

change_of_basis(b::PeriodicTranslationDict, ::Type{LinearCombinationOfPeriodicTranslationDict}; options...) = DualGram(b; options...)

function coefficients_in_other_basis{B<:LinearCombinationOfPeriodicTranslationDict}(b::PeriodicTranslationDict, ::Type{B}; options...)
    e = zeros(b)
    e[1] = 1
    change_of_basis(b, B; options...)*e
end

BasisFunctions.extension_operator(s1::LinearCombinationOfPeriodicTranslationDict, s2::LinearCombinationOfPeriodicTranslationDict; options...) =
    wrap_operator(s1, s2, change_of_basis(s2; options...)*extension_operator(superdict(s1), superdict(s2))*inv(change_of_basis(s1; options...)))

BasisFunctions.restriction_operator(s1::LinearCombinationOfPeriodicTranslationDict, s2::LinearCombinationOfPeriodicTranslationDict; options...) =
    wrap_operator(s1, s2, change_of_basis(s2; options...)*restriction_operator(superdict(s1), superdict(s2))*inv(change_of_basis(s1; options...)))

"""
  Set representing the dual basis.
"""
struct DualPeriodicTranslationDict{T} <: LinearCombinationOfPeriodicTranslationDict{PeriodicTranslationDict, T}
    superdict       :: PeriodicTranslationDict{T}
    coefficients    :: Array{T,1}
end



DualPeriodicTranslationDict(set::PeriodicTranslationDict{T}; options...) where {T} =
    DualPeriodicTranslationDict{T}(set, coefficients_in_other_basis(set, LinearCombinationOfPeriodicTranslationDict; options...))

BasisFunctions.superdict(b::DualPeriodicTranslationDict) = b.superdict

BasisFunctions.dual(b::DualPeriodicTranslationDict; options...) = superdict(b)

BasisFunctions.Gram(b::DualPeriodicTranslationDict; options...) = inv(Gram(superdict(b); options...))

"""
  Set representing the dual basis with respect to a discrete norm on the oversampled grid.
"""
struct DiscreteDualPeriodicTranslationDict{T} <: LinearCombinationOfPeriodicTranslationDict{PeriodicTranslationDict, T}
    superdict       :: PeriodicTranslationDict{T}
    coefficients    :: Array{T,1}

    oversampling    :: T
end



function DiscreteDualPeriodicTranslationDict(set::PeriodicTranslationDict{T}; oversampling=BasisFunctions.default_oversampling(set), options...) where {T}
    DiscreteDualPeriodicTranslationDict{T}(set, coefficients_in_other_basis(set, DiscreteDualPeriodicTranslationDict; oversampling=oversampling, options...), oversampling)
end

BasisFunctions.superdict(b::DiscreteDualPeriodicTranslationDict) = b.superdict
BasisFunctions.coefficients(b::DiscreteDualPeriodicTranslationDict) = b.coefficients

BasisFunctions.default_oversampling(b::DiscreteDualPeriodicTranslationDict) = b.oversampling

BasisFunctions.dual(b::DiscreteDualPeriodicTranslationDict; options...) = superdict(b)

change_of_basis(b::PeriodicTranslationDict, ::Type{DiscreteDualPeriodicTranslationDict}; options...) = DiscreteDualGram(b; options...)

BasisFunctions.resize(b::DiscreteDualPeriodicTranslationDict, n::Int) = discrete_dual(resize(dual(b), n); oversampling=BasisFunctions.default_oversampling(b))

BasisFunctions.dict_promote_domaintype(b::DiscreteDualPeriodicTranslationDict{T}, ::Type{S}) where {T,S} =
    discrete_dual(promote_domaintype(dual(b), S); oversampling=BasisFunctions.default_oversampling(b))
