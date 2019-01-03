##################
# Platform
##################

# 1D generators
primal_bspline_generator(::Type{T}, degree::Int) where {T} = n->BSplineTranslatesBasis(n, degree, T)
struct DualBSplineGenerator
    primal_generator
    oversampling::Int
    DualBSplineGenerator(primal_generator, oversampling::Int) = (@assert BasisFunctions.isdyadic(oversampling); new(primal_generator, oversampling))
end

function (DG::DualBSplineGenerator)(n::Int)
    B = DG.primal_generator(n)
    DG = DiscreteDualGram(B, oversampling=DG.oversampling)
    OperatedDict(DG)
end

dual_bspline_generator(primal_bspline_generator, oversampling::Int) = DualBSplineGenerator(primal_bspline_generator, oversampling)

# ND generators
primal_bspline_generator(::Type{T}, degree1::Int, degree2::Int, degree::Int...) where {T} = primal_bspline_generator(T, [degree1, degree2, degree...])

primal_bspline_generator(::Type{T}, degree::NTuple{N,Int}) where {N,T} = tensor_generator(T, map(d->primal_bspline_generator(T, d), degree)...)

(DG::DualBSplineGenerator)(n1::Int, n2::Int, ns::Int...) = DG([n1,ns...])

function (DG::DualBSplineGenerator)(n::Array{Int})
    B = DG.primal_generator(n)
    DG = DiscreteDualGram(B, oversampling=DG.oversampling)
    tensorproduct([OperatedDict(DGi) for DGi in elements(DG)]...)
end

# Sampler
bspline_sampler(::Type{T}, primal, oversampling::Int) where {T} = n-> GridSamplingOperator(gridbasis(grid(primal(n*oversampling)),T))

# params
bspline_param(init::Int) = DoublingSequence(init)

bspline_param(init::NTuple{N,Int}) where N = TensorSequence([MultiplySequence(i,2 .^(1/length(init))) for i in init])

# Platform
bspline_platform(::Type{T}, init::NTuple{1,Int}, degree::NTuple{1,Int}, oversampling::Int) where {T} =
    bspline_platform(T, init[1], degree[1], oversampling)
function bspline_platform(::Type{T}, init::Union{Int,NTuple{N,Int}}, degree::Union{Int,NTuple{N,Int}}, oversampling::Int) where {N,T}
	primal = primal_bspline_generator(T, degree)
	dual = dual_bspline_generator(primal, oversampling)
        sampler = bspline_sampler(T, primal, oversampling)
        dual_sampler = n->(1/length(dual(n)))*sampler(n)
	params = bspline_param(init)
	GenericPlatform(primal = primal, dual = dual, sampler = sampler, dual_sampler=dual_sampler,
		params = params, name = "B-Spline translates")
end

function bspline_platform(::Type{T}, init::Int, degree::Int, oversampling::Int) where {T}
	primal = primal_bspline_generator(T, degree)
	dual = dual_bspline_generator(primal, oversampling)
        sampler = bspline_sampler(T, primal, oversampling)
        dual_sampler = n->(1/length(dual(n)))*sampler(n)
	params = bspline_param(init)
	GenericPlatform(primal = primal, dual = dual, sampler = sampler, dual_sampler=dual_sampler,
		params = params, name = "B-Spline translates")
end
