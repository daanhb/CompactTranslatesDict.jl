support(b::CompactTranslationDict, idx) = support(b, native_index(b, idx))

support(b::CompactTranslationDict{T}, idx::TransIndex) where {T} = periodize_interval(value(idx)*step(b)+kernel_support(b), support(b))

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
        UnionDomain((Interval(infimum(outdomain), supremum(indomain)-per)),(Interval(infimum(indomain), supremum(outdomain))))
    end
end

if !isdefined(BasisFunctions,:default_eval_expansion)
    function default_eval_expansion(dict::Dictionary, coefficients, x)
        @assert size(coefficients) == size(dict)

        T = span_codomaintype(dict)
        z = zero(T)
        # It is safer below to use eval_element than unsafe_eval_element, because of
        # the check on the support.
        @inbounds for idx in eachindex(coefficients)
            z = z + coefficients[idx] * eval_element(dict, idx, x)
        end
        z
    end
else
    using BasisFunctions: default_eval_expansion
end


compactlysupported(D::FullSpace) = false
compactlysupported(D::AbstractHyperSphere) = true
compactlysupported(D::DomainSets.MappedDomain) = compactlysupported(superdomain(D))
compactlysupported(D::AbstractInterval) = true
compactlysupported(D::ProductDomain) = all(map(compactlysupported, elements(D))...)
