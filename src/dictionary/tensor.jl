
const CompactTranslatesTensorProductDict{N,TUPLE,S,T} = BasisFunctions.TensorProductDict{N,TUPLE,S,T} where {N,TUPLE<: Tuple{Vararg{CompactTranslatesDict.CompactTranslationDict}},S,T}
const BSplineTensorProductDict{N,TUPLE,S,T} = BasisFunctions.TensorProductDict{N,TUPLE,S,T} where {N,TUPLE<: Tuple{Vararg{CompactTranslatesDict.DiffPeriodicBSplineBasis}},S,T}

degree(T::BSplineTensorProductDict) = map(degree, elements(T))
Bdiff(T::BSplineTensorProductDict) = map(Bdiff, elements(T))

grid_evaluation_operator(s::CompactTranslatesTensorProductDict, dgs::GridBasis, grid::ProductGrid; options...) =
    tensorproduct([grid_evaluation_operator(si, dgsi, gi; options...) for (si, dgsi, gi) in zip(elements(s), elements(dgs), elements(grid))]...)
