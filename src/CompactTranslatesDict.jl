module CompactTranslatesDict

using BasisFunctions, DomainSets, FillArrays,
    GridArrays, InfiniteVectors, Reexport

include("util/PeriodicIntervals.jl")
@reexport using .PeriodicIntervals

include("util/infinitevectorscompat.jl")
include("util/compactinfinitevector.jl")

include("translates.jl")
include("PeriodicEquispacedTranslatesDicts.jl")
include("DiffPeriodicBSplineBases.jl")
include("CompactPeriodicEquispacedTranslatesDuals.jl")

end
