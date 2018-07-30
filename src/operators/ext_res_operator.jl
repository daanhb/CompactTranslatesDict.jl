const My1DIndexType = Union{Int,AbstractVector{Int}}
if VERSION < v"0.7-"
    const MyNDIndexType{N} = Union{CartesianIndex{N},AbstractVector{CartesianIndex{N}},CartesianRange{CartesianIndex{N}}}
else
    const MyNDIndexType{N} = Union{CartesianIndex{N},AbstractVector{CartesianIndex{N}},CartesianIndices{N}}
end
const MyIndexType = Union{My1DIndexType,MyNDIndexType}

struct ExtResOperator{ELT} <: BasisFunctions.DictionaryOperator{ELT}
    src::Dictionary
    dest::Dictionary
    srcindices
    destindices
    op::DictionaryOperator{ELT}
end

function ExtResOperator(destindices::My1DIndexType, op::DictionaryOperator, srcindices::My1DIndexType)
    @assert is_fastly_indexable(op)
    ExtResOperator(src(op)[srcindices], dest(op)[destindices], srcindices, destindices, op)
end

function ExtResOperator(destindices::MyNDIndexType, op::TensorProductOperator, srcindices::MyNDIndexType)
    @assert VERSION < v"0.7-" ? reduce(&, true, map(is_fastly_indexable, elements(op))) : reduce(&, map(is_fastly_indexable, elements(op)); init=true) 
    ExtResOperator(src(op)[srcindices], dest(op)[destindices], srcindices, destindices, op)
end

ExtResOperator(op::DictionaryOperator, srcindices::MyIndexType) =
    ExtResOperator(eachindex(dest(op)), op, srcindices)

ExtResOperator(destindices::MyIndexType, op::DictionaryOperator) =
    ExtResOperator(destindices, op, eachindex(src(op)))

# ExtResOperator(res::IndexRestrictionOperator, op::DictionaryOperator, ext::IndexExtensionOperator) =
#     ExtResOperator(subindices(res), op, subindices(ext))
#
# ExtResOperator(op::DictionaryOperator, ext::IndexExtensionOperator) =
#     ExtResOperator(op, subindices(ext))
#
# ExtResOperator(res::IndexRestrictionOperator, op::DictionaryOperator) =
#     ExtResOperator(subindices(res), op)


getindex(M::DictionaryOperator, i::My1DIndexType, j::My1DIndexType) =
    ExtResOperator(i, M, j)

getindex(M::DictionaryOperator, i::Colon, j::My1DIndexType) =
    ExtResOperator(M, j)

getindex(M::DictionaryOperator, i::My1DIndexType, j::Colon) =
    ExtResOperator(i, M)

getindex(M::ExtResOperator, i::My1DIndexType, j::My1DIndexType) =
    ExtResOperator(M.destindices[i], M.op, M.srcindices[j])

getindex(M::ExtResOperator, i::Colon, j::My1DIndexType) =
    ExtResOperator(M.destindices, M.op, M.srcindices[j])

getindex(M::ExtResOperator, i::My1DIndexType, j::Colon) =
    ExtResOperator(M.destindices[i], M.op, M.srcindices)

getindex(M::TensorProductOperator, i::MyNDIndexType, j::MyNDIndexType) =
    ExtResOperator(i, M, j)

getindex(M::TensorProductOperator, i::Colon, j::MyNDIndexType) =
    ExtResOperator(M, j)

getindex(M::TensorProductOperator, i::MyNDIndexType, j::Colon) =
    ExtResOperator(i, M)

function getindex(M::TensorProductOperator, i::My1DIndexType, j::My1DIndexType)
    warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    l = native_index.(src(M), j)
    ExtResOperator(k, M, l)
end

function getindex(M::TensorProductOperator, i::Colon, j::My1DIndexType)
    warn("Deprecated: partial linear indexing of multidimensional object")
    l = native_index.(src(M), j)
    ExtResOperator(M, l)
end

function getindex(M::TensorProductOperator, i::My1DIndexType, j::Colon)
    warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    ExtResOperator(k, M)
end



getindex(M::ExtResOperator, i::Int, j::Int) =
    fast_getindex(M.op, M.destindices[i], M.srcindices[j])

matrix(M::ExtResOperator) =
    fast_matrix(M, M.op, M.srcindices, M.destindices)

fast_matrix(M::ExtResOperator, T::TensorProductOperator, srcindices::MyNDIndexType, destindices::MyNDIndexType) =
    fast_matrix!(_zeros(M), elements(T), srcindices, destindices)

_zeros(M::ExtResOperator) =
    _zeros(eltype(M), size(M), M.destindices, M.srcindices)

_zeros(::Type{T}, size, destindices, srcindices) where {T} =
    zeros(T, size)

_zeros(::Type{T}, size, destindices::Union{CartesianIndex,Int}, srcindices) where {T} =
    zeros(T, length(srcindices))

_zeros(::Type{T}, size, destindices, srcindices::Union{CartesianIndex,Int}) where {T} =
    zeros(T, length(destindices))

fast_matrix(M::ExtResOperator, T::DictionaryOperator, srcindices::My1DIndexType, destindices::My1DIndexType) =
    [fast_getindex(T, i, j) for i in destindices, j in srcindices]

function fast_matrix!(m::Matrix, ops, srcindices, destindices)
    for (i,k) in enumerate(destindices), (j, l) in enumerate(srcindices)
        m[i,j] = fast_getindex(ops, k, l)
    end
    m
end

function fast_matrix!(m::Vector, ops, srcindices, destindices::Union{CartesianIndex,Int})
    for (i,l) in enumerate(srcindices)
        m[i] = fast_getindex(ops, destindices, l)
    end
    m
end

function fast_matrix!(m::Vector, ops, srcindices::Union{CartesianIndex,Int}, destindices)
    for (i,k) in enumerate(destindices)
        m[i] = fast_getindex(ops, k, srcindices)
    end
    m
end

fast_getindex(M::TensorProductOperator, i::CartesianIndex, j::CartesianIndex) =
    fast_getindex(elements(M), i, j)

@generated function fast_getindex(ops, i::CartesianIndex{N}, j::CartesianIndex{N}) where {N}
    l = [:(fast_getindex(ops[$i], i[$i], j[$i])) for i in 1:N]
    # println("compile fast_getindex for tensor dimension $N")
    s = ""
    for i in 1:N-1
        s *= string(l[i])*"*"
    end
    s *= string(l[end])
    return parse(s)
end
