const My1DIndexType = Union{AbstractVector{Int}}
const NdOperator{ELT} = Union{TensorProductOperator{ELT},OperatorSum{ELT}}
const MyNDIndexType{N} = Union{AbstractVector{CartesianIndex{N}},CartesianIndices{N}}
const MyIndexType = Union{My1DIndexType,MyNDIndexType}

struct ExtResOperator{ELT} <: BasisFunctions.DictionaryOperator{ELT}
    src::Dictionary
    dest::Dictionary
    srcextension    ::IndexExtensionOperator
    destrestriction ::IndexRestrictionOperator
    op::DictionaryOperator
    RopE::CompositeOperator
    function ExtResOperator{ELT}(src, dest, srcindices, destindices, op) where ELT
        @assert !isa(src, Subdictionary)
        @assert !isa(dest, Subdictionary)
        if size(srcindices) != size(src)
            res_src = src[srcindices]
        else
            res_src = src
        end
        if size(destindices) != size(dest)
            res_dest = dest[destindices]
        else
            res_dest = dest
        end
        E = IndexExtensionOperator(res_src, src, srcindices)
        R = IndexRestrictionOperator(dest, res_dest, destindices)
        new{ELT}(res_src, res_dest, E, R, op, R*op*E)
    end
end

srcindices(M::ExtResOperator) = subindices(M.srcextension)
destindices(M::ExtResOperator) = subindices(M.destrestriction)

function ExtResOperator(destindices::My1DIndexType, op::DictionaryOperator{ELT}, srcindices::My1DIndexType) where{ELT}
    @assert is_fastly_indexable(op)
    ExtResOperator{ELT}(src(op), dest(op), srcindices, destindices, op)
end

function ExtResOperator(destindices::MyNDIndexType, op::NdOperator{ELT}, srcindices::MyNDIndexType) where {ELT}
    @assert is_fastly_indexable(op)
    ExtResOperator{ELT}(src(op), dest(op), srcindices, destindices, op)
end

apply!(M::ExtResOperator, a, b) = apply!(M.RopE, a, b)

unsafe_wrap_operator(src, dest, op::ExtResOperator{ELT}) where ELT =
    ExtResOperator{ELT}(full_dict(src), full_dict(dest), srcindices(op),  destindices(op), op.op)

full_dict(d::Subdictionary) = superdict(d)
full_dict(d::Dictionary) = d
full_dict(d::GridBasis) = gridbasis(full_grid(grid(d)), coefficienttype(d))
full_grid(a::AbstractGrid) = a
full_grid(a::IndexSubGrid) = supergrid(a)

adjoint(M::ExtResOperator{ELT}) where ELT = ExtResOperator{ELT}(full_dict(dest(M)), full_dict(src(M)), destindices(M), srcindices(M), M.op')

is_fastly_indexable(::BasisFunctions.ArrayOperator) = true

fast_getindex(m::BasisFunctions.ArrayOperator, i::Int, j::Int) = getindex(m.A,i,j)

is_fastly_indexable(m::ExtResOperator) = is_fastly_indexable(m.A)
is_fastly_indexable(m::AbstractArray) = true

is_fastly_indexable(a...) = false

is_fastly_indexable(T::TensorProductOperator) =
    reduce(&, map(is_fastly_indexable, elements(T)); init=true)
is_fastly_indexable(T::OperatorSum) =
    reduce(&, map(is_fastly_indexable, elements(T)); init=true)

tovector(a::AbstractArray) = reshape(a, length(a))
ExtResOperator(op::DictionaryOperator, srcindices::MyIndexType) =
    ExtResOperator(tovector(eachindex(dest(op))), op, srcindices)

ExtResOperator(destindices::MyIndexType, op::DictionaryOperator) =
    ExtResOperator(destindices, op, tovector(eachindex(src(op))))

ExtResOperator(op::DictionaryOperator) =
    ExtResOperator(tovector(eachindex(dest(op))), op, tovector(eachindex(src(op))))

getindex(M::DictionaryOperator, i::My1DIndexType, j::My1DIndexType) =
    ExtResOperator(i, M, j)

getindex(M::DictionaryOperator, i::Colon, j::My1DIndexType) =
    ExtResOperator(M, j)

getindex(M::DictionaryOperator, i::My1DIndexType, j::Colon) =
    ExtResOperator(i, M)

getindex(M::DictionaryOperator, i::Colon, j::Colon) =
    ExtResOperator(M)

getindex(M::ExtResOperator, i::My1DIndexType, j::My1DIndexType) =
    ExtResOperator(destindices(M)[i], M.op, srcindices(M)[j])

getindex(M::ExtResOperator, i::Colon, j::My1DIndexType) =
    ExtResOperator(destindices(M), M.op, srcindices(M)[j])

getindex(M::ExtResOperator, i::My1DIndexType, j::Colon) =
    ExtResOperator(destindices(M)[i], M.op, srcindices(M))

# TODO implement more of these
function getindex(M::ExtResOperator, i::MyNDIndexType, j::My1DIndexType)
    @assert length(size(destindices(M))) > 1
    ExtResOperator(destindices(M)[i], M.op, srcindices(M)[j])
end

function getindex(M::ExtResOperator, i::My1DIndexType, j::MyNDIndexType)
    @assert length(size(srcindices(M))) > 1
    ExtResOperator(destindices(M)[i], M.op, srcindices(M)[i])
end

getindex(M::NdOperator, i::MyNDIndexType, j::MyNDIndexType) =
    ExtResOperator(i, M, j)

getindex(M::NdOperator, i::Colon, j::MyNDIndexType) =
    ExtResOperator(M, j)

getindex(M::NdOperator, i::MyNDIndexType, j::Colon) =
    ExtResOperator(i, M)

function getindex(M::NdOperator, i::My1DIndexType, j::My1DIndexType)
    @warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    l = native_index.(src(M), j)
    ExtResOperator(k, M, l)
end

function getindex(M::NdOperator, i::Colon, j::My1DIndexType)
    @warn("Deprecated: partial linear indexing of multidimensional object")
    l = native_index.(src(M), j)
    ExtResOperator(M, l)
end

function getindex(M::NdOperator, i::My1DIndexType, j::Colon)
    @warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    ExtResOperator(k, M)
end


# checkbounds(op::DictionaryOperator, i::Colon, j::Int) =
# 	1 <= j <= size(op,2) || throw(BoundsError())
#
# checkbounds(op::DictionaryOperator, i::Int, j::Colon) =
# 	1 <= i <= size(op,1) || throw(BoundsError())
# checkbounds(op::DictionaryOperator, i::Colon, j::Colon) = true

unsafe_getindex(M::ExtResOperator, i::Int, j::Int) =
    fast_getindex(M.op, destindices(M)[i], srcindices(M)[j])

matrix(M::ExtResOperator) =
    fast_matrix(M, M.op, srcindices(M), destindices(M))

matrix!(M::ExtResOperator, m::Array) =
    fast_matrix!(m, M.op, srcindices(M), destindices(M))

_zeros(M::ExtResOperator) =
    _zeros(eltype(M), size(M), destindices(M), srcindices(M))

_zeros(::Type{T}, size, destindices, srcindices) where {T} =
    zeros(T, size)

_zeros(::Type{T}, size, destindices::Union{CartesianIndex,Int}, srcindices) where {T} =
    zeros(T, length(srcindices))

_zeros(::Type{T}, size, destindices, srcindices::Union{CartesianIndex,Int}) where {T} =
    zeros(T, length(destindices))

fast_matrix(M::ExtResOperator, T::DictionaryOperator, srcindices::My1DIndexType, destindices::My1DIndexType) =
    [fast_getindex(T, i, j) for i in destindices, j in srcindices]

fast_matrix(M::ExtResOperator, T::NdOperator, srcindices::MyNDIndexType, destindices::MyNDIndexType) =
    fast_matrix!(_zeros(M), T, srcindices, destindices)

function fast_matrix!(m::AbstractMatrix, T, srcindices, destindices)
    for (i,k) in enumerate(destindices), (j, l) in enumerate(srcindices)
        m[i,j] = fast_getindex(T, k, l)
    end
    m
end

fast_getindex(T, k::CartesianIndex{1}, l::CartesianIndex{1}) = fast_getindex(T, k[1], l[1])

fast_matrix!(m::AbstractMatrix, T::NdOperator, srcindices, destindices) =
    fast_matrix!(m, elements(T), T, srcindices, destindices)

@noinline function fast_matrix!(m::AbstractMatrix, ops, T::TensorProductOperator, srcindices, destindices)
    for (i,k) in enumerate(destindices), (j, l) in enumerate(srcindices)
        m[i,j] = fast_getindex(ops, T, k, l)
    end
    m
end

fast_matrix!(m::AbstractMatrix, Sops, S::OperatorSum, srcindices, destindices) =
    fast_matrix!(m, map(elements, Sops), Sops, coefficients(S), S, srcindices, destindices)

@noinline function fast_matrix!(m::AbstractMatrix, TSops, Sops, coefficients, S::OperatorSum, srcindices, destindices)
    for (i,k) in enumerate(destindices), (j, l) in enumerate(srcindices)
        m[i,j] = fast_getindex(TSops, Sops, coefficients, S, k, l)
    end
    m
end

function fast_matrix!(m::AbstractVector, T::NdOperator, srcindices, destindices::Union{CartesianIndex,Int})
    for (i,l) in enumerate(srcindices)
        m[i] = fast_getindex(T, destindices, l)
    end
    m
end

function fast_matrix!(m::AbstractVector, T::NdOperator, srcindices::Union{CartesianIndex,Int}, destindices)
    for (i,k) in enumerate(destindices)
        m[i] = fast_getindex(T, k, srcindices)
    end
    m
end

@generated function fast_getindex(ops, T::TensorProductOperator, i::CartesianIndex{N}, j::CartesianIndex{N}) where {N}
    l = [:(fast_getindex(ops[$i], i[$i], j[$i])) for i in 1:N]
    # println("compile fast_getindex for tensor dimension $N")
    s = ""
    for i in 1:N-1
        s *= string(l[i])*"*"
    end
    s *= string(l[end])
    return Meta.parse(s)
end

@generated function fast_getindex(TSops, Sops::Tuple{Vararg{DictionaryOperator,M}}, coefficients, op::OperatorSum, i::CartesianIndex{N}, j::CartesianIndex{N}) where {N,M}
    l = [:(coefficients[$i]*fast_getindex(TSops[$i], Sops[$i], i, j)) for i in 1:M]
    s = ""
    for i in 1:M-1
        s *= string(l[i])*"+"
    end
    s *= string(l[end])
    return Meta.parse(s)
end
