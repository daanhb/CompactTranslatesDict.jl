
export compactinfinitevector
"""
    compactinfinitevector(dict::Dictionary, grid::AbstractGrid)

Return a CompactInfiniteVector `vec` that holds the evaluation of dict in grid.
We assume that dict has the structure of a PeriodicEquispacedTranslates dictionary, i.e, the
evaluation matrix is a VerticalBandedOperator using an AbstractEquispacedGrid grid.

In multiple dimensions it returns a tuple of vectors: (vec_1, ..., vec_D)

# Example
```jldocs
julia> dict = BSplineTranslatesBasis(10,3);
julia> g = interpolation_grid(dict);
julia> v = compactinfinitevector(dict, g);

julia> vv = v[1:length(g)];
julia> vv .+= v[-length(g)+1:0];
julia> @test vv ≈ evaluation_matrix(dict[1], g)
Test Passed
```
"""
function compactinfinitevector(dict::Dictionary{T}, grid::AbstractEquispacedGrid; verbose=false, threshold=sqrt(eps(real(T))), options...) where {T}
    @assert rem(length(grid),length(dict)) == 0
    A = evaluation(dict, grid; verbose=verbose, threshold=threshold, options...)
    @assert A isa VerticalBandedOperator
    a = convert(Vector{T}, A.A.array)
    f = 1
    l = length(a)
    for i in 1:length(a)
        if !isapprox(a[i] + 1, 1, rtol=threshold, atol=threshold)
            break
        end
        f += 1
    end
    for j in length(a):-1:1
        if !isapprox(a[j] + 1, 1, rtol=threshold, atol=threshold)
            break
        end
        l -= 1
    end

    truncatedarray = a[f:l]
    os = mod(A.A.offset+f, length(grid))
    if !(0 <= os + length(truncatedarray) <= length(grid))
        os -= length(grid)
    end
    CompactInfiniteVector(truncatedarray, os)
end
