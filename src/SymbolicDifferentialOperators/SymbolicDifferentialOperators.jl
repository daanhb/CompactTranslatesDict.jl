module SymbolicDifferentialOperators

export order

import Base.string

alphabet = "abcdefghijklmnopqrstuvwxyz"

abstract type AbstractDiffOperator end
order(diffoperator::AbstractDiffOperator) = diffoperator.order

dictionary(::AbstractDiffOperator) = Dict{Char,Int}()

import Base.show
show(io::IO, d::AbstractDiffOperator) = show(io, string(d))

abstract type AbstractDiffCoefficient end

abstract type AbstractScalarDiffCoefficient <: AbstractDiffCoefficient end

scalar(s::AbstractScalarDiffCoefficient) = s.scalar

abstract type AbstractDiffOperator1d{DN}  <: AbstractDiffOperator where {DN} end

dimension_name(diffoperator::AbstractDiffOperator1d{DN}) where DN  = DN

dimension_names(d::AbstractDiffOperator1d) = (dimension_name(d),)

struct IdentityCoefficient <: AbstractScalarDiffCoefficient
end

scalar(::IdentityCoefficient) = 1

import Base.==

==(c::Number, ::IdentityCoefficient) = (c==1)
==(::IdentityCoefficient, c::Number) = (c==1)

struct ConstantCoefficient{T<:Number} <: AbstractScalarDiffCoefficient
    scalar::T
end

string(c::ConstantCoefficient) = string(scalar(c))

==(c::Number, C::ConstantCoefficient) = (c==scalar(C))
==(C::ConstantCoefficient, c::Number) = (c==scalar(C))
==(C::ConstantCoefficient, I::IdentityCoefficient) = (scalar(C)==1)
==(I::IdentityCoefficient, C::ConstantCoefficient) = (scalar(C)==1)

string(c::IdentityCoefficient) = "1"

struct IdentityDifferentialOperator <: AbstractDiffOperator
end

I = IdentityDifferentialOperator()

assigned_number(::IdentityDifferentialOperator) = 0

order(::IdentityDifferentialOperator) = 0

string(::IdentityDifferentialOperator) = "1"

slice(::IdentityDifferentialOperator, x::Char, dns::NTuple{N,Char}) where N = IdentityDifferentialOperator()

is_first_dimension_name(d::AbstractDiffOperator, x::Char, dns::NTuple{N,Char}=dimension_names(d)) where N =
    x == dns[1]

struct ZeroDifferentialOperator <: AbstractDiffOperator
end

O = ZeroDifferentialOperator()

assigned_number(::ZeroDifferentialOperator) = -1

order(::ZeroDifferentialOperator) = 0

string(::ZeroDifferentialOperator) = "0"

slice(::ZeroDifferentialOperator, x::Char, dns::NTuple{N,Char}) where N = ZeroDifferentialOperator()

struct ScaledDifferentialOperator{C<:AbstractDiffCoefficient,D<:AbstractDiffOperator} <: AbstractDiffOperator
    coeff::C
    diff::D
end

assigned_number(s::ScaledDifferentialOperator) = assigned_number(operator(s))

coefficient(s::ScaledDifferentialOperator) = s.coeff

operator(s::ScaledDifferentialOperator) = s.diff

dimension_names(s::ScaledDifferentialOperator) = dimension_names(operator(s))

slice(s::ScaledDifferentialOperator, x::Char, dns::NTuple{N,Char}=dimension_names(s)) where N = is_first_dimension_name(operator(s), x, dns) ?
    coefficient(s)*slice(operator(s), x, dns) :
    1*slice(operator(s), x, dns)


==(s1::ScaledDifferentialOperator, s2::ScaledDifferentialOperator) =
    (coefficient(s1) == coefficient(s2) && (operator(s1) == operator(s2)))
==(d::AbstractDiffOperator, s::ScaledDifferentialOperator) =
    (coefficient(s) == 1 && (operator(s) == d))
==(s::ScaledDifferentialOperator, d::AbstractDiffOperator) =
    (coefficient(s) == 1 && (operator(s) == d))

string(s::ScaledDifferentialOperator) = string(coefficient(s))*string(operator(s))
string(s::ScaledDifferentialOperator{C,IdentityDifferentialOperator})  where C<:AbstractDiffCoefficient = string(coefficient(s))
string(s::ScaledDifferentialOperator{IdentityCoefficient,IdentityDifferentialOperator}) = "1"
string(s::ScaledDifferentialOperator{IdentityCoefficient,D}) where D<:AbstractDiffOperator= string(operator(s))

ScaledandUnscaledDifferentialOperators{D} = Union{D, ScaledDifferentialOperator{C, D} where  C<:AbstractDiffCoefficient} where {D<:AbstractDiffOperator}
SaUDIffOpp{D} = ScaledandUnscaledDifferentialOperators{D} where D<:AbstractDiffOperator

struct PartialDifferentialOperator{DN} <: AbstractDiffOperator1d{DN}
    order::Int
end

for c in alphabet
    delta = Symbol("δ_"*string(c))
    D = Symbol("D_"*string(c))
    deltas = Symbol("δ"*string(c))
    Ds = Symbol("D"*string(c))
    @eval begin
        $delta = PartialDifferentialOperator{$c}(1)
        $D = PartialDifferentialOperator{$c}(1)
        $deltas = PartialDifferentialOperator{$c}(1)
        $Ds = PartialDifferentialOperator{$c}(1)
        export $delta, $D, $deltas, $Ds
    end
end

dictionary(p::PartialDifferentialOperator) = Dict(dimension_name(p)=>order(p))

is_first_dimension_name(p::PartialDifferentialOperator, x::Char, ::NTuple{N,Char}) where N =
    (x==dimension_name(p))

slice(p::PartialDifferentialOperator, x::Char, ::NTuple{N,Char}=dimension_names(p)) where N = (x==dimension_name(p)) ?
    p :
    ZeroDifferentialOperator()

import Base:^
^(p::PartialDifferentialOperator, o::Int) = PartialDifferentialOperator{dimension_name(p)}(order(p)*o)

struct ProductDifferentialOperator{DNS} <: AbstractDiffOperator
    orders::NTuple{N,Int} where N
    function ProductDifferentialOperator{DNS}(values) where DNS
        @assert length(DNS) == length(values)
        p = sortperm(collect(DNS))
        new{DNS[p]}(values[p])
    end
end

==(p1::Union{PartialDifferentialOperator,ProductDifferentialOperator}, p2::Union{PartialDifferentialOperator,ProductDifferentialOperator}) =
    dictionary(p1) == dictionary(p2)

order(p::ProductDifferentialOperator) = p.orders

dimension_names(p::ProductDifferentialOperator{DNS}) where DNS = DNS

dictionary(p::ProductDifferentialOperator) = Dict(k=>v for (k,v) in zip(dimension_names(p), order(p)))

slice(p::ProductDifferentialOperator, x::Char, ::NTuple{N,Char}=dimension_names(p)) where N = (x∈dimension_names(p)) ?
     PartialDifferentialOperator{x}(order(p)[find(x.==dimension_names(p))[1]]) :
     ZeroDifferentialOperator()

is_first_dimension_name(p::ProductDifferentialOperator, x::Char, ::NTuple{N,Char}) where N =
    (x==dimension_names(p)[1])

DiffOps = Union{SaUDIffOpp{ProductDifferentialOperator{DNS}} where DNS, SaUDIffOpp{PartialDifferentialOperator{DN}} where DN}

function ProductDifferentialOperator(partialoperators::DiffOps...)
    c = 1
    d = Dict{Char,Int}()
    for p in partialoperators
        if isa(p, ScaledDifferentialOperator)
            c *=  coefficient(p)
            p = operator(p)
        end
        if isa(p, PartialDifferentialOperator)
            add_order!(d, dimension_name(p), order(p))
        else
            for (k,v) in zip(dimension_names(p), order(p))
                add_order!(d, k, v)
            end
        end
    end
    c*ProductDifferentialOperator{tuple(keys(d)...)}(tuple(values(d)...))
end

function string(p::Union{PartialDifferentialOperator,ProductDifferentialOperator})
    r = ""
    for (k,v) in dictionary(p)
        r *= "δ"*string(k)*string(superscript(v))
    end
    r
end

function add_order!(d::Dict, key::Char, value::Int)
    if haskey(d, key)
        old_order = pop!(d, key)
        push!(d, key=>old_order+value)
    else
        push!(d, key=>value)
    end
end

function superscript(i::Int)
    if i > 9
        join(superscript(d) for d in reverse(digits(i)))
    elseif i==1
        Char(185)
    elseif i==2
        Char(178)
    elseif i==3
        Char(179)
    elseif i < 0
        error()
    else
        Char(8304+i)
    end
end

function assigned_number(p::Union{PartialDifferentialOperator,ProductDifferentialOperator})
    dns = dimension_names(p)
    a = 0
    for (c,o) in dictionary(p)
        for i in 1:o
            a += Int(c)-Int('a')+1
            a = a << 5
        end
    end
    a = a >> 5
end

import Base: *
*(c::Number, d::AbstractDiffOperator) = ConstantCoefficient(c)*d
*(c::AbstractDiffCoefficient, d::AbstractDiffOperator) = (c==1) ?
    ScaledDifferentialOperator(IdentityCoefficient(), d) :
    ScaledDifferentialOperator(c, d)
*(diffops::DiffOps...) = ProductDifferentialOperator(diffops...)
*(s::Number, c::AbstractDiffCoefficient) = ConstantCoefficient(s)*c
*(s::IdentityCoefficient, ::IdentityCoefficient) = IdentityCoefficient()
*(::IdentityCoefficient, c::ConstantCoefficient) = c
*(c::ConstantCoefficient, ::IdentityCoefficient) = c
*(c1::ConstantCoefficient, c2::ConstantCoefficient) = ConstantCoefficient(scalar(c1)*scalar(c2))
*(d::AbstractDiffOperator, s::ScaledDifferentialOperator) = coefficient(s)*(d*operator(s))
*(s::ScaledDifferentialOperator, d::AbstractDiffOperator) = coefficient(s)*(operator(s)*d)
*(s1::ScaledDifferentialOperator, s2::ScaledDifferentialOperator) =
    (coefficient(s1)*coefficient(s2))*(operator(s1)*operator(s2))
*(c::AbstractDiffCoefficient, s::ScaledDifferentialOperator) = (c*coefficient(s))*operator(s)
*(s::ScaledDifferentialOperator, c::AbstractDiffCoefficient) = (c*coefficient(s))*operator(s)
*(d::AbstractDiffOperator, ::IdentityDifferentialOperator) = d
*(::IdentityDifferentialOperator, d::AbstractDiffOperator) = d

struct SumDifferentialOperator{DNS} <: AbstractDiffOperator
    ops::Vector{ScaledDifferentialOperator}
    order::NTuple{N,Int} where N

    function SumDifferentialOperator{DNS}(ops, values) where DNS
        @assert length(DNS) == length(values)
        p = sortperm(collect(DNS))
        new{DNS[p]}(ops,values[p])
    end
end

function SumDifferentialOperator(ops::ScaledDifferentialOperator...)
    V = collect(ops)
    A = ScaledDifferentialOperator[]
    v = assigned_number.(V)
    t = nothing
    d = Dict{Char,Int}()
    for i in unique(v)
        c = 0
        for vi in V
            if assigned_number(vi) == i
                t = operator(vi)
                c +=  coefficient(vi)
            end
        end
        if !isa(t, ZeroDifferentialOperator)
            for (k,v) in dictionary(t)
                max_order!(d, k, v)
            end
            push!(A, c*t)
        end
    end
    sort!(A, by=assigned_number)
    if length(A) > 1
        SumDifferentialOperator{tuple(keys(d)...)}(A, tuple(values(d)...))
    else
        A[1]
    end
end


function max_order!(d::Dict, key::Char, value::Int)
    if haskey(d, key)
        old_value = pop!(d, key)
        push!(d, key=>max(old_value,value))
    else
        push!(d, key=>value)
    end
end

elements(s::SumDifferentialOperator) = s.ops

dimension_names(::SumDifferentialOperator{DNS}) where DNS = DNS

slice(s::SumDifferentialOperator, x::Char, dns::NTuple{N,Char}=dimension_names(s)) where N =
    SumDifferentialOperator([slice(e, x, dns) for e in elements(s)]...)

function string(p::SumDifferentialOperator)
    r = ""
    r *= string(elements(p)[1])
    for i in 2:length(elements(p))
        r *= "+"
        r *= string(elements(p)[i])
    end
    r
end

==(s1::SumDifferentialOperator, s2::SumDifferentialOperator) =
    elements(s1) == elements(s2)

import Base.+
+(c::Number, d::AbstractDiffOperator) = c*I+d
+(d::AbstractDiffOperator, c::Number) = c*I+d
+(c1::AbstractScalarDiffCoefficient, c2::AbstractDiffCoefficient) = ConstantCoefficient(scalar(c1)+scalar(c2))
+(c1::AbstractScalarDiffCoefficient, c2::Number) = c1 + ConstantCoefficient(c2)
+(c1::Number, c2::AbstractScalarDiffCoefficient) = ConstantCoefficient(c1) + c2
+(d1::AbstractDiffOperator, d2::AbstractDiffOperator) =
    d1 == d2 ?
        ConstantCoefficient(2)*d1 :
        SumDifferentialOperator(1*d1, 1*d2)
+(s1::ScaledDifferentialOperator, s2::ScaledDifferentialOperator) =
    operator(s1) == operator(s2) ?
        (coefficient(s1)+coefficient(s2))*operator(s1) :
        SumDifferentialOperator(s1, s2)
+(s1::SumDifferentialOperator, s2::AbstractDiffOperator) =
    s1 + 1*s2
+(s1::AbstractDiffOperator, s2::SumDifferentialOperator) =
    1*s1 + s2
+(s1::SumDifferentialOperator, s2::ScaledDifferentialOperator) =
    SumDifferentialOperator(elements(s1)..., s2)
+(s1::ScaledDifferentialOperator, s2::SumDifferentialOperator) =
    SumDifferentialOperator(s1, elements(s2)...)
+(s1::SumDifferentialOperator, s2::SumDifferentialOperator) =
    SumDifferentialOperator(elements(s1)..., elements(s2)...)
+(::ZeroDifferentialOperator, d::AbstractDiffOperator) = d
+(d::AbstractDiffOperator, ::ZeroDifferentialOperator) = d
+(d::ZeroDifferentialOperator, ::ZeroDifferentialOperator) = d
+(::ScaledDifferentialOperator{C,ZeroDifferentialOperator} where C<:AbstractDiffCoefficient, d::AbstractDiffOperator) = d
+(d::AbstractDiffOperator, ::ScaledDifferentialOperator{C,ZeroDifferentialOperator} where C<:AbstractDiffCoefficient) = d
+(::ScaledDifferentialOperator{C,ZeroDifferentialOperator} where C<:AbstractDiffCoefficient, ::ScaledDifferentialOperator{C,ZeroDifferentialOperator} where C<:AbstractDiffCoefficient) = ZeroDifferentialOperator()
import Base.-
-(c::Number, d2::AbstractDiffOperator) = d2 + (-c)
-(d2::AbstractDiffOperator, c::Number) = d2 + (-c)
-(d1::AbstractDiffOperator, d2::AbstractDiffOperator) = d1 + ConstantCoefficient(-1)*d2

*(c::AbstractDiffCoefficient, s::SumDifferentialOperator) =
    SumDifferentialOperator([c*e for e in elements(s)]...)
*(c::AbstractDiffOperator, s::SumDifferentialOperator) =
    SumDifferentialOperator([c*e for e in elements(s)]...)
*(c::ScaledDifferentialOperator, s::SumDifferentialOperator) =
    SumDifferentialOperator([c*e for e in elements(s)]...)
*(s1::SumDifferentialOperator,s2::SumDifferentialOperator) =
    SumDifferentialOperator([e1*e2 for e1 in elements(s1) for e2 in elements(s2)]...)
end
