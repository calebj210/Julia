using LinearAlgebra
using AssociatedLegendrePolynomials

function sphericalharmonic(p::NamedTuple, m::I, l::I) where I <: Integer
    if m == 0
        sph = Plm(l, m, cos(p.θ))
    elseif m < 0
        m = -m
        coeff = (-1)^m * sqrt(2factorial(big(l - m)) / factorial(big(l + m)))
        sph = coeff * Plm(l, m, cos(p.θ)) * sin(m * p.ϕ)
    else
        coeff = (-1)^m * sqrt(2factorial(big(l - m)) / factorial(big(l + m)))
        sph = coeff * Plm(l, m, cos(p.θ)) * cos(m * p.ϕ)
    end

    return p.r^l * sph
end

function sphericalharmonic(x::T, y::T, z::T, d::I, n::I) where {T <: Number, I <: Integer}
    p = (;
        r = norm([x,y,z]),
        ϕ = atan(y, x),
        θ = atan(norm([x,y]), z)
    )

    return sphericalharmonic(p, d, n)
end

sphericalharmonic(x::Union{Vector,Tuple}, d, n) = sphericalharmonic(x..., d, n)

function sphericalharmonic(x, N::Integer)
    deg = floor(Int64, sqrt(N) + 1) - 1
    N -= (deg)^2
    m = -deg + N

    return sphericalharmonic(x, m, deg)
end

function circularharmonic(p::NamedTuple, N::Integer)
    d = N ÷ 2
    if N == 0
        return 1
    elseif N % 2 == 1
        return p.r^(d + 1) * cos((d + 1) * p.θ)
    else
        return p.r^d * sin(d * p.θ)
    end

end

circularharmonic(x::T, y::T, N::Integer) where T <: Number = circularharmonic((; r = norm([x,y]), θ = atan(y, x)), N)
circularharmonic(x::Union{Vector,Tuple}, N::Integer) = circularharmonic(x..., N)

harmonic(node::NTuple{3,<:Number}, N::Integer) = sphericalharmonic(node, N)
harmonic(node::NTuple{2,<:Number}, N::Integer) = circularharmonic(node, N)
