if VERSION < v"0.7-"
    Pkg.clone("https://github.com/vincentcp/CardinalBSplines.git")
    Pkg.clone("https://github.com/daanhb/BasisFunctions.jl.git")
    Pkg.checkout("BasisFunctions", "julia-0.7")
    Pkg.build("BasisFunctions")
    Pkg.clone("https://github.com/daanhb/Domains.jl.git")
end
