function generate_bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::SumDifferentialOperator, dns::NTuple{N,Char}=dimension_names(d)) where N where T
    dicts = tuple(map(x->generate_bspline_basis(T, size, degree, operator(x), dns), SymbolicDifferentialEquations.elements(d))...)
    coefficients = [T(scalar(coefficient(x))) for x in SymbolicDifferentialEquations.elements(d)]
    CompactTranslatesDict.CompactTranslationDictSum(dicts, coefficients)
end

function generate_bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::PartialDifferentialOperator, dns::NTuple{N,Char}) where {N,T}
    diff = zeros(Int,length(size))
    diff[find(dimension_name(d).==dns)] = order(d)
    TensorProductDict([DiffBSplineTranslatesBasis(size[i], degree[i], diff[i], T) for i in 1:length(size)]...)
end

function generate_bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::ProductDifferentialOperator, dns::NTuple{N,Char}) where {N,T}
    diff = zeros(Int,length(size))
    for (i,c) in enumerate(dns)
        for (x,v) in SymbolicDifferentialEquations.dictionary(d)
            if c == x
                diff[i] = v
            end
        end
    end
    TensorProductDict([DiffBSplineTranslatesBasis(size[i], degree[i], diff[i], T) for i in 1:length(size)]...)
end

function generate_bspline_basis(::Type{T}, size::NTuple{N,Int}, degree::NTuple{N,Int}, d::IdentityDifferentialOperator, dns::NTuple{N,Char}) where {N,T}
    TensorProductDict([DiffBSplineTranslatesBasis(size[i], degree[i], 0, T) for i in 1:length(size)]...)
end
