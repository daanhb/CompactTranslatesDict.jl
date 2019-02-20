
struct BSplinePlatform{T,D} <: BasisPlatform
end

BSplinePlatform(degree::Int=3) = BSplinePlatform{Float64,degree}()

dictionary(p::BSplinePlatform{T,D}, n) where {T,D} = BSplineTranslatesBasis{T,D}(n)

first_parameters(p::BSplinePlatform) = (8,8)

SolverStyle(p::BSplinePlatform, ::OversamplingStyle) = DualStyle()

measure(p::BSplinePlatform{T}) where {T} = FourierMeasure{T}()



struct BSplineExtensionPlatform{T,D} <: FramePlatform
    basisplatform   ::  BSplinePlatform{T,D}
    domain          ::  Domain

    function BSplineExtensionPlatform(basisplatform::BSplinePlatform{T,D}, domain::Domain{T}) where {T,D}
        @assert issubset(domain, support(measure(basisplatform)))
        new{T,D}(basisplatform, domain)
    end
end

BSplineExtensionPlatform(domain::Domain{T}, D::Int=3) where {T} =
    BSplineExtensionPlatform(BSplinePlatform{T,D}(), domain)

dictionary(p::BSplineExtensionPlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

measure(platform::BSplineExtensionPlatform) = measure(dictionary(platform, 1))
