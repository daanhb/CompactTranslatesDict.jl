

struct CompactTranslationDictSum{N,S,T} <: DerivedDict{S,T}
    superdict::CompactTranslatesTensorProductDict{N,DT,S,T} where DT
    dicts::Tuple{Vararg{Dictionary,M}} where M
    coefficients::Vector{T}

    function CompactTranslationDictSum{N,S,T}(dicts, coefs) where {N,S,T}
        superdict = dicts[1]
        if VERSION < v"0.7-"
            @assert reduce(&, true, size(superdict)==size(dict) for dict in dicts)
            @assert reduce(&, true, infimum(support(superdict))≈infimum(support(dict)) for dict in dicts)
            @assert reduce(&, true, supremum(support(superdict))≈supremum(support(dict)) for dict in dicts)
        else
            @assert reduce(&, size(superdict)==size(dict) for dict in dicts; init=true)
            @assert reduce(&, infimum(support(superdict))≈infimum(support(dict)) for dict in dicts; init=true)
            @assert reduce(&, supremum(support(superdict))≈supremum(support(dict)) for dict in dicts; init=true)
        end
        @assert eltype(coefs) == T
        @assert length(dicts) == length(coefs)
        new{N,S,T}(superdict, tuple(dicts...), coefs)
    end
end

CompactTranslationDictSum(dicts::Tuple{Vararg{CompactTranslatesTensorProductDict{N,TUPLE,S,T} where {TUPLE}}}, coefficients::Vector{T} = ones(T, length(dicts))) where {N,S,T} =
    CompactTranslationDictSum{N,S,T}(dicts, coefficients)

elements(S::CompactTranslationDictSum) = S.dicts

element(S::CompactTranslationDictSum, i) = S.dicts[i]

coefficients(S::CompactTranslationDictSum) = S.coefficients

function BasisFunctions.promote_coefficienttype(S::CompactTranslationDictSum, ::Type{T}) where T
    new_dicts = map(promote_coefficienttype, elements(S))
    new_coef_type = promote(eltype((coefficients(S))), T)
    new_coefs = zeros(T, size(coefficients((S))))
    new_coefs .= coefficients((S))
    CompactTranslatesDict(new_dicts, new_coefs)
end

# function similar_dictionary(s::CompactTranslationDictSum, s2::Dictionary)
#     @show s
#     @show s2
#     CompactTranslationDictSum(s2)
# end

unsafe_eval_element(dict::CompactTranslationDictSum, i::ProductIndex, x) =
    sum(c*unsafe_eval_element(e, i, x) for (c,e) in zip(coefficients(dict),elements(dict)))

oversampled_grid(dict::CompactTranslationDictSum, sampling_factor::Real) = ProductGrid([oversampled_grid(e, sampling_factor) for e in elements(dict.superdict)]...)

resize(dict::CompactTranslationDictSum, s) = CompactTranslationDictSum(tuple([resize(e, s) for e in elements(dict)]...), coefficients(dict))

grid_evaluation_operator(S::CompactTranslationDictSum, dgs::GridBasis, grid::ProductGrid; options...) =
    OperatorSum(S, dgs, tuple([tensorproduct([grid_evaluation_operator(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)
        for s in elements(S)]...), coefficients(S))

grid_evaluation_operator(s::CompactTranslationDictSum, dgs::GridBasis, subgrid::AbstractSubGrid; options...) =
    BasisFunctions._grid_evaluation_operator(s, dgs, subgrid; options...)

# function UnNormalizedGram(s::CompactTranslationDictSum, oversampling = 1)
#     grid = oversampled_grid(s, oversampling)
#     TensorCirculantOperator(evaluation_operator(s, grid)'*evaluation_operator(s, grid))
# end
