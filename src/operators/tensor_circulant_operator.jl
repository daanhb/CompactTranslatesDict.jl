# tensor_circulant_operator.jl

"""
A tensor circulant operator is a operator circulant in every dimension.

Several other operators can be converted into a circulant matrix, and this
conversion happens automatically when such operators are combined into a composite
operator.
"""
struct TensorCirculantOperator{T,N} <: DictionaryOperator{T} where N
    src                 :: Dictionary
    dest                :: Dictionary
    F                   :: DictionaryOperator
    iF                  ::DictionaryOperator
    eigenvaluematrix    :: AbstractArray{C,N} where C<:Complex
    Fscratch
    iFscratch

    TensorCirculantOperator{T,N}(op_src::Dictionary, op_dest::Dictionary, F::DictionaryOperator, iF::DictionaryOperator, matrix::AbstractArray) where {T,N} =
      new{T,N}(op_src, op_dest, F, iF, matrix, zeros(dest(F)), zeros(dest(iF)))
end

src(c::TensorCirculantOperator) = c.src
dest(c::TensorCirculantOperator) = c.dest

TensorCirculantOperator(src::Dictionary, dest::Dictionary, D::AbstractArray; options...) = TensorCirculantOperator(complex(eltype(D)), src, dest, fft(D); options...)

function TensorCirculantOperator(::Type{T}, op_src::Dictionary, op_dest::Dictionary, opD::AbstractArray{T,N}; real_circulant_tol=sqrt(eps(real(T))), verbose=false, options...) where {T,N}
    cpx_src = TensorProductDict([DiscreteVectorDictionary{complex(eltype(opD))}(size(opD, i)) for i in 1:N]...)
    A = promote_type(eltype(opD),T)

    F = TensorProductOperator([forward_fourier_operator(element(cpx_src, i), element(cpx_src, i), A; verbose=verbose, options...) for i in 1:N]...)
    iF = inv(F)
    #realify a circulant operator if src and dest are real (one should imply the other).
    if isreal(op_src) && isreal(op_dest)
        imag_norm = norm(imag(fft(opD)))
        imag_norm > real_circulant_tol && warn("realified circulant operator, lost an accuracy of $(imag_norm)")
        r_S, r_D, r_A = op_eltypes(op_src, op_dest, real(T))
        r_src = promote_coeftype(op_src, r_S)
        r_dest = promote_coeftype(op_dest, r_D)

        return TensorCirculantOperator{r_A,N}(r_src, r_dest, F, iF, opD)

    end
    TensorCirculantOperator{A,N}(op_src, op_dest, F, iF, opD)
end

function TensorCirculantOperator(op::DictionaryOperator{T}) where {T}
    e = zeros(T, src(op))
    e[1] = one(T)
    C = TensorCirculantOperator(src(op), dest(op), op*e)
    e = map(T,rand(src(op)))

    @assert C*e≈op*e
    C
end

inner_array(C::TensorCirculantOperator) = C.eigenvaluematrix

similar_operator(op::TensorCirculantOperator, src, dest) =
    TensorCirculantOperator(src, dest, op.eigenvaluematrix)

inv(C::TensorCirculantOperator{A,N}) where {A,N} = TensorCirculantOperator{A,N}(dest(C), src(C), C.F, C.iF, 1 ./ C.eigenvaluematrix)

function element_wise_pinv(a::AbstractArray, tol)
    r = pinv.(a)
    r[abs.(r) .< tol] .= 0
    r
end

# What tolerance should be used for the pinv here?
pinv(C::TensorCirculantOperator{T}, tolerance = eps(numtype(C))) where {T} = TensorCirculantOperator(src(c), dest(c), element_wise_pinv(C.eigenvaluematrix, tolerance))

for op in (:+, :-, :*)
    @eval $op(c1::TensorCirculantOperator{A,N}, c2::TensorCirculantOperator{A,N}) where {A,N} = TensorCirculantOperator{A,N}(src(c2), dest(c1), c1.F, c2.iF, broadcast($op,inner_array(c1),inner_array(c2)))
end

*(scalar::Real, c::TensorCirculantOperator{A,N}) where {A,N} = TensorCirculantOperator{A,N}(src(c), dest(c), c.F, c.iF, scalar*inner_array(c))
*(c::TensorCirculantOperator, scalar::Real) = scalar*c

apply!(c::TensorCirculantOperator, coef_dest, coef_src) = apply!(c, coef_dest, coef_src, Val{isreal(c)})
function apply!(c::TensorCirculantOperator, coef_dest, coef_src, ::Type{Val{false}})
    apply!(c.F, c.Fscratch, coef_src)
    c.Fscratch .*= c.eigenvaluematrix
    apply!(c.iF, coef_dest, c.Fscratch)
end

function apply!(c::TensorCirculantOperator, coef_dest, coef_src, ::Type{Val{true}})
    apply!(c.F, c.Fscratch, coef_src)
    c.Fscratch .*= c.eigenvaluematrix
    apply!(c.iF, c.iFscratch, c.Fscratch)
    copy!(coef_dest, real(c.iFscratch))
end


apply_inplace!(c::TensorCirculantOperator, coef_dest) = apply_inplace!(c, coef_dest, Val{isreal(c)})
function apply_inplace!(c::TensorCirculantOperator, coef_dest, ::Type{Val{false}})
    copy!(c.Fscratch, coef_dest)
    apply_inplace!(c.F, c.Fscratch)
    c.Fscratch .*= c.eigenvaluematrix
    apply_inplace!(c.iF, c.Fscratch)
    copy!(coef_dest, c.Fscratch)
end

function apply_inplace!(c::TensorCirculantOperator, coef_dest, ::Type{Val{true}})
    copy!(c.Fscratch, coef_dest)
    apply_inplace!(c.F, c.Fscratch)
    c.Fscratch ./= c.eigenvaluematrix
    apply_inplace!(c.iF, c.Fscratch)
    copy!(coef_dest, c.Fscratch)
end

# TODO find to do this without copying code (TensorCirculantOperator is no DerivedOperator?)
function matrix!(op::TensorCirculantOperator, a)
    coef_src  = zeros(src(op))
    coef_dest = zeros(dest(op))
    matrix_fill!(op, a, coef_src, coef_dest)
end

# TODO find to do this without copying code (TensorCirculantOperator is no DerivedOperator?)
function unsafe_getindex(op::TensorCirculantOperator, i::Int, j::Int)
  coef_src = zeros(src(op))
	coef_dest = zeros(dest(op))
	coef_src[j] = one(eltype(op))
	apply!(op, coef_dest, coef_src)
	coef_dest[i]
end
