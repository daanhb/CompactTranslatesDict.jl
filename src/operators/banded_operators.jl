# banded_oprators.jl

"""
A banded operator of which every row contains equal elements.

The top row starts at index offset, the second row at step+offset.
"""
struct IndexableHorizontalBandedOperator{ELT} <: BasisFunctions.DictionaryOperator{ELT}
    src::Dictionary
    dest::Dictionary
    array::Vector{ELT}
    step::Int
    offset::Int
    function IndexableHorizontalBandedOperator{ELT}(src::Dictionary, dest::Dictionary, array::Vector{ELT}, step::Int=1, offset::Int=0) where ELT
        @assert length(array) <= length(src)
        @assert step <= length(src) # apply! only works if step is smaller then L
        new{ELT}(src, dest, array, step, offset)
    end
end

IndexableHorizontalBandedOperator(src::Dictionary, dest::Dictionary, array::Vector{ELT}, step::Int=1, offset::Int=0) where ELT =
    IndexableHorizontalBandedOperator{ELT}(src, dest, array, step, offset)

fast_getindex(op::IndexableHorizontalBandedOperator, i::Int, j::Int) =
    _get_horizontal_banded_index(op.array, op.step, op.offset, size(op,1), size(op,2), i, j)

function _get_horizontal_banded_index(array::Vector{ELT}, step::Int, offset::Int, M::Int, N::Int, i::Int, j::Int) where {ELT}
    # first transform to an index of the first column
    index = mod(j-step*(i-1)-1-offset, N)+1
    if index <= length(array)
        array[index]
    else
        ELT(0)
    end
end

matrix(M::IndexableHorizontalBandedOperator) = [getindex(M, i, j) for i in 1:size(M,1), j in 1:size(M,2)]

function apply!(op::IndexableHorizontalBandedOperator, dest::Vector, src::Vector)
    dest[:] .= 0
    L = length(src)
    aL = length(op.array)
    dL = length(dest)
    @inbounds for a_i in 1:aL
        ind = mod(a_i+op.offset-1,L)+1
        for d_i in 1:dL
            dest[d_i] += op.array[a_i]*src[ind]
            ind += op.step
            if ind > L
                ind -= L
            end
        end
    end
    dest
end

is_fastly_indexable(::IndexableHorizontalBandedOperator) = true
getindex(M::IndexableHorizontalBandedOperator, i::Int, j::Int) = fast_getindex(M, i, j)

adjoint(M::IndexableHorizontalBandedOperator) = IndexableVerticalBandedOperator(dest(M), src(M), M.array, M.step, M.offset)

struct IndexableVerticalBandedOperator{ELT} <: BasisFunctions.DictionaryOperator{ELT}
    src::Dictionary
    dest::Dictionary
    array::Vector{ELT}
    step::Int
    offset::Int
    function IndexableVerticalBandedOperator{ELT}(src::Dictionary, dest::Dictionary, array::Vector{ELT}, step::Int, offset::Int) where ELT
        @assert length(array) <= length(dest)
        @assert step <= length(dest) # apply! only works if step is smaller then L
        new{ELT}(src, dest, array, step, offset)
    end
end

IndexableVerticalBandedOperator(src::Dictionary, dest::Dictionary, array::Vector{ELT}, step::Int=1, offset::Int=0) where ELT =
    IndexableVerticalBandedOperator{ELT}(src, dest, array, step, offset)

fast_getindex(op::IndexableVerticalBandedOperator, i::Int, j::Int) =
    _get_vertical_banded_index(op.array, op.step, op.offset, size(op,1), size(op,2), i, j)

function apply!(op::IndexableVerticalBandedOperator, dest::Vector, src::Vector)
    dest[:] .= 0
    # assumes step is smaller then L
    L = length(dest)
    aL = length(op.array)
    sL = length(src)
    @inbounds for a_i in 1:aL
        ind = mod(a_i+op.offset-1,L)+1
        for s_i in 1:sL
            dest[ind] += op.array[a_i]*src[s_i]
            ind += op.step
            if ind > L
                ind -= L
            end
        end
    end

    dest
end

function _get_vertical_banded_index(array::Vector{ELT}, step::Int, offset::Int, M::Int, N::Int, i::Int, j::Int) where {ELT}
    # first transform to an index of the first column
    index = mod(i-step*(j-1)-1-offset, M)+1
    if index <= length(array)
        array[index]
    else
        ELT(0)
    end
end

matrix(M::IndexableVerticalBandedOperator) = [getindex(M, i, j) for i in 1:size(M,1), j in 1:size(M,2)]

is_fastly_indexable(::IndexableVerticalBandedOperator) = true
getindex(M::IndexableVerticalBandedOperator, i::Int, j::Int) = fast_getindex(M, i, j)

adjoint(M::IndexableVerticalBandedOperator) = IndexableHorizontalBandedOperator(dest(M), src(M), M.array, M.step, M.offset)
