module CompactTranslatesDict

using Reexport
include("PeriodicIntervals.jl")
@reexport using .PeriodicIntervals

include("InfiniteVectorsCompat.jl")
using .InfiniteVectorsCompat

include("TranslatesDictionaries.jl/TranslatesDictionaries.jl")
@reexport using .TranslatesDictionaries

include("TranslatesDictionaries.jl/PeriodicEquispacedTranslatesDicts.jl")
@reexport using .PeriodicEquispacedTranslatesDicts

include("TranslatesDictionaries.jl/DiffPeriodicBSplineBases.jl")
@reexport using .DiffPeriodicBSplineBases

include("CompactInfiniteVectors.jl")
@reexport using .CompactInfiniteVectors

include("TranslatesDictionaries.jl/CompactPeriodicEquispacedTranslatesDuals.jl")
@reexport using .CompactPeriodicEquispacedTranslatesDuals

end
