function bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::SumDifferentialOperator, dns::NTuple{N,Char}=dimension_names(d)) where N where T
    dicts = tuple(map(x->bspline_basis(T, size, degree, operator(x), dns), SymbolicDifferentialOperators.elements(d))...)
    coefficients = [T(scalar(coefficient(x))) for x in SymbolicDifferentialOperators.elements(d)]
    CompactTranslatesDict.CompactTranslationDictSum(dicts, coefficients)
end

function bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::PartialDifferentialOperator, dns::NTuple{N,Char}) where {N,T}
    diff = zeros(Int,length(size))
    diff[find(dimension_name(d).==dns)] = order(d)
    tensor_product_dict(size, degree, tuple(diff...), T)
end

function bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::ProductDifferentialOperator, dns::NTuple{N,Char}) where {N,T}
    diff = zeros(Int,length(size))
    for (i,c) in enumerate(dns)
        for (x,v) in SymbolicDifferentialOperators.dictionary(d)
            if c == x
                diff[i] = v
            end
        end
    end
    tensor_product_dict(size, degree, tuple(diff...), T)
end

bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::IdentityDifferentialOperator, dns::NTuple{N,Char}) where {N,T} =
    tensor_product_dict(size, degree, tuple(diff...), T)

bspline_evaluation_operator(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, grid::ProductGrid, d::SymbolicDifferentialOperators.AbstractDiffOperator, dns::NTuple{N,Char}=dimension_names(d)) where {N,T} =
    evaluation_operator(bspline_basis(T, size, degree, d, dns), grid)
