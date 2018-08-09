
"A linear combination of operators."
struct OperatorSum{T,S,N} <: DictionaryOperator{T}
    src::Dictionary
    dest::Dictionary
    ops ::Tuple{Vararg{DictionaryOperator,N}}
    coefficients::  Vector{T}
    scratch     ::  S

    function OperatorSum{T,S,N}(s, d, ops, coefficients, scratch) where {T,S,N}
        # We don't enforce that source and destination of op1 and op2 are the same, but at least
        # their sizes and coefficient types must match.
        if VERSION < v"0.7-"
            @assert reduce(&, true, size(src(ops[1]))==size(src(op)) for op in ops)
            @assert reduce(&, true, size(dest(ops[1]))==size(dest(op)) for op in ops)
            @assert reduce(&, true, T==coeftype(dest(op)) for op in ops)
            @assert reduce(&, true, T==coeftype(dest(op)) for op in ops)
        else 
            @assert reduce(&, size(src(ops[1]))==size(src(op)) for op in ops; init=true)
            @assert reduce(&, size(dest(ops[1]))==size(dest(op)) for op in ops; init=true)
            @assert reduce(&, T==coeftype(dest(op)) for op in ops; init=true)
            @assert reduce(&, T==coeftype(dest(op)) for op in ops; init=true)
        end
        @assert length(ops) == length(coefficients)
        new{T,S,N}(s, d, ops, coefficients, scratch)
    end
end

# function OperatorSum(ops, coefficients)
#     s = CompactTranslationDictSum(map(src, ops), coefficients)
#     d = CompactTranslationDictSum(map(dest, ops), coefficients)
#     OperatorSum(s, d, ops, coefficients)
# end

function OperatorSum(s, d, ops, coefficients)
    T = promote_type(map(eltype, ops)...)
    scratch = zeros(dest(ops[1]))
    OperatorSum{T,typeof(scratch),length(ops)}(s, d, ops, [T(c) for c in coefficients], scratch)
end

elements(op::OperatorSum) = op.ops

element(op::OperatorSum, i::Int) = elements(op)[i]

coefficients(op::OperatorSum) = op.coefficients

function apply!(op::OperatorSum{S,T,N}, dest, src, coef_dest, coef_src) where {S,T,N}
    coef_dest[:] .= 0

    for j in 1:N
        apply!(element(op, j), op.scratch, coef_src)
        c = op.coefficients[j]
        for i in eachindex(coef_dest)
            coef_dest[i] += c * op.scratch[i]
        end
    end
    coef_dest
end

adjoint(op::OperatorSum) = OperatorSum(dest(op), src(op), adjoint.(elements(op)), conj.(coefficients(op)))
