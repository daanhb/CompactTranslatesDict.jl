const My1DIndexType = Union{AbstractVector{Int}}
const NdOperator{ELT} = Union{TensorProductOperator{ELT},OperatorSum{ELT}}

if VERSION < v"0.7-"
    const MyNDIndexType{N} = Union{AbstractVector{CartesianIndex{N}},CartesianRange{CartesianIndex{N}}}
else
    const MyNDIndexType{N} = Union{AbstractVector{CartesianIndex{N}},CartesianIndices{N}}
end
const MyIndexType = Union{My1DIndexType,MyNDIndexType}

struct ExtResOperator{ELT} <: BasisFunctions.DictionaryOperator{ELT}
    src::Dictionary
    dest::Dictionary
    srcextension    ::IndexExtensionOperator
    destrestriction ::IndexRestrictionOperator
    op::DictionaryOperator{ELT}
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
        new(res_src, res_dest, E, R, op, R*op*E)
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


is_fastly_indexable(::ExtResOperator) = true

is_fastly_indexable(a...) = false

is_fastly_indexable(T::TensorProductOperator) = VERSION < v"0.7-" ?
    reduce(&, true, map(is_fastly_indexable, elements(T))) :
    reduce(&, map(is_fastly_indexable, elements(T)); init=true)
is_fastly_indexable(T::OperatorSum) = VERSION < v"0.7-" ?
    reduce(&, true, map(is_fastly_indexable, elements(T))) :
    reduce(&, map(is_fastly_indexable, elements(T)); init=true)

ExtResOperator(op::DictionaryOperator, srcindices::MyIndexType) =
    ExtResOperator(eachindex(dest(op)), op, srcindices)

ExtResOperator(destindices::MyIndexType, op::DictionaryOperator) =
    ExtResOperator(destindices, op, eachindex(src(op)))

ExtResOperator(op::DictionaryOperator) =
    ExtResOperator(eachindex(dest(op)), op, eachindex(src(op)))

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

if VERSION < v"0.7-"
    getindex(cr::CartesianRange, i::AbstractVector{Int}) = collect(cr)[i]
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
    warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    l = native_index.(src(M), j)
    ExtResOperator(k, M, l)
end

function getindex(M::NdOperator, i::Colon, j::My1DIndexType)
    warn("Deprecated: partial linear indexing of multidimensional object")
    l = native_index.(src(M), j)
    ExtResOperator(M, l)
end

function getindex(M::NdOperator, i::My1DIndexType, j::Colon)
    warn("Deprecated: partial linear indexing of multidimensional object")
    k = native_index.(dest(M), i)
    ExtResOperator(k, M)
end



getindex(M::ExtResOperator, i::Int, j::Int) =
    fast_getindex(M.op, destindices(M)[i], srcindices(M)[j])

matrix(M::ExtResOperator) =
    fast_matrix(M, M.op, srcindices(M), destindices(M))

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

function fast_matrix!(m::Matrix, T::NdOperator, srcindices, destindices)
    for (i,k) in enumerate(destindices), (j, l) in enumerate(srcindices)
        m[i,j] = fast_getindex(T, k, l)
    end
    m
end

function fast_matrix!(m::Vector, T::NdOperator, srcindices, destindices::Union{CartesianIndex,Int})
    for (i,l) in enumerate(srcindices)
        m[i] = fast_getindex(T, destindices, l)
    end
    m
end

function fast_matrix!(m::Vector, T::NdOperator, srcindices::Union{CartesianIndex,Int}, destindices)
    for (i,k) in enumerate(destindices)
        m[i] = fast_getindex(T, k, srcindices)
    end
    m
end

fast_getindex(M::TensorProductOperator, i::CartesianIndex, j::CartesianIndex) =
    tensor_fast_getindex(elements(M), i, j)

@generated function tensor_fast_getindex(ops, i::CartesianIndex{N}, j::CartesianIndex{N}) where {N}
    l = [:(fast_getindex(ops[$i], i[$i], j[$i])) for i in 1:N]
    # println("compile fast_getindex for tensor dimension $N")
    s = ""
    for i in 1:N-1
        s *= string(l[i])*"*"
    end
    s *= string(l[end])
    return parse(s)
end

fast_getindex(M::OperatorSum, i::CartesianIndex, j::CartesianIndex) =
    sum_fast_getindex(elements(M), coefficients(M), i, j)

@generated function sum_fast_getindex(ops::Tuple{Vararg{DictionaryOperator,M}}, coefficients, i::CartesianIndex{N}, j::CartesianIndex{N}) where {N,M}
    l = [:(coefficients[$i]*fast_getindex(ops[$i], i, j)) for i in 1:M]
    s = ""
    for i in 1:M-1
        s *= string(l[i])*"+"
    end
    s *= string(l[end])
    return parse(s)
end
