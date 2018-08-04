

struct CompactTranslationDictSum{N,S,T} <: DerivedDict{S,T}
    superdict::CompactTranslatesTensorProductDict{N,DT,S,T} where DT
    dicts::Tuple{Vararg{Dictionary,M}} where M
    coefficients::Vector{T}

    function CompactTranslationDictSum{N,S,T}(dicts, coefs) where {N,S,T}
        superdict = dicts[1]
        @assert reduce(&, true, size(superdict)==size(dict) for dict in dicts)
        @assert reduce(&, true, infimum(support(superdict))≈infimum(support(dict)) for dict in dicts)
        @assert reduce(&, true, supremum(support(superdict))≈supremum(support(dict)) for dict in dicts)
        @assert eltype(coefs) == T
        @assert length(dicts) == length(coefs)
        new{N,S,T}(superdict, tuple(dicts...), coefs)
    end
end

CompactTranslationDictSum(dicts::Tuple{Vararg{CompactTranslatesTensorProductDict{N,TUPLE,S,T} where {TUPLE}}}, coefficients::Vector{T} = ones(T, length(dicts))) where {N,S,T} =
    CompactTranslationDictSum{N,S,T}(dicts, coefficients)

elements(S::CompactTranslationDictSum) = S.dicts

coefficients(S::CompactTranslationDictSum) = S.coefficients

unsafe_eval_element(dict::CompactTranslationDictSum, i::ProductIndex, x) =
    sum(c*unsafe_eval_element(e, i, x) for (c,e) in zip(coefficients(dict),elements(dict)))

resize(dict::CompactTranslationDictSum, s) = (@assert(size(dict)==s); dict)

grid_evaluation_operator(S::CompactTranslationDictSum, dgs::GridBasis, grid::ProductGrid; options...) =
    OperatorSum(tuple([tensorproduct([grid_evaluation_operator(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)
                for s in elements(S)]...), coefficients(S))
