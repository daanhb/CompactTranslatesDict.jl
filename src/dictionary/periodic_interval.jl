
"""
For two given intervals `[a,b]` and `[A,B]`, a periodic interval represents
the intersection of `[A,B]` with the periodic repetition of `[a,b]` with period
`B-A`:

`[a +k*(B-A),b+k*(B-A))] ∩ [A,B]`.

Depending on the relative location of the two intervals, the periodic interval
may be a single interval or a union of two intervals.
"""
struct PeriodicInterval{I1,I2,T} <: Domain{T}
    subdomain       ::  I1
    periodicdomain  ::  I2
end

PeriodicInterval(subdomain::AbstractInterval{T}, periodicdomain::AbstractInterval{T}) where {T} =
    PeriodicInterval{typeof(subdomain),typeof(periodicdomain),T}(subdomain, periodicdomain)

DomainSets.indomain(x, d::PeriodicInterval) = _indomain(x, d, d.subdomain, d.periodicdomain)

function _indomain(x, d::PeriodicInterval, subdomain::AbstractInterval, periodicdomain::AbstractInterval)
    if x ∉ periodicdomain
        # x is not contained in [A,B]
        return false
    end

    a,b = extrema(subdomain)
    A,B = extrema(periodicdomain)
    if b-a >= B-A
        # interval [a,b] is larger than [A,B]
        return true
    end
    if a >= A && b <= B
        # interval [a,b] is contained within [A,B]
        return x ∈ subdomain
    else
        # wrap [a,b] to the union of [a1,B] and [A,b1]
        a1 = A + mod(a-A,B-A)
        b1 = A + mod(b-A,B-A)
        return (x >= a1) || (x <= b1)
    end
end

approx_indomain(x, d::PeriodicInterval, tolerance) = _approx_indomain(x, d, tolerance, d.subdomain, d.periodicdomain)

function _approx_indomain(x, d::PeriodicInterval, tolerance, subdomain::AbstractInterval, periodicdomain::AbstractInterval)
    if !approx_in(x, periodicdomain)
        # x is not contained in [A,B]
        return false
    end

    a,b = extrema(subdomain)
    A,B = extrema(periodicdomain)
    if b-a >= B-A
        # interval [a,b] is larger than [A,B]
        return true
    end
    if a >= A && b <= B
        # interval [a,b] is contained within [A,B]
        return approx_in(x, subdomain)
    else
        # wrap [a,b] to the union of [a1,B] and [A,b1]
        a1 = A + mod(a-A,B-A)
        b1 = A + mod(b-A,B-A)
        return (x >= a1-tolerance) || (x <= b1+tolerance)
    end
end

function infimum(d::PeriodicInterval)
    a, b = extrema(d.subdomain)
    A, B = extrema(d.periodicdomain)
    a1 = A + mod(a-A,B-A)
    b1 = A + mod(b-A,B-A)
    if b1 > a1
        a1
    else
        A
    end
end

function extremum(d::PeriodicInterval)
    a, b = extrema(d.subdomain)
    A, B = extrema(d.periodicdomain)
    a1 = A + mod(a-A,B-A)
    b1 = A + mod(b-A,B-A)
    if b1 > a1
        b1
    else
        B
    end
end
