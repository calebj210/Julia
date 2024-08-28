#=
#   ODE approach to computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: August 28, 2024
=#

using Polynomials
using SpecialFunctions
include("TimeStep.jl")

function sgn(x)
    iszero(x) ? -one(x) : sign(x)
end

function pfq_taylor(a, b, z, N)
    coeffs = zeros(typeof(z), N + 1)
    coeffs[1] = 1.0

    for n ∈ 1 : N
        coeffs[n + 1] = coeffs[n] * prod(a .+ (n - 1)) / prod(b .+ (n - 1)) / n 
    end

    p  = Polynomial(coeffs)
    dp = derivative(p)

    return [p(z), dp(z)]
end

function F21(a, b, c, z::Number; h = .1, z0 = 0im, order = 20, taylorN = 100)
    if iszero(z0)
        z0 = sgn(imag(z)) * 0.5im
    end

    y0, dy0 = pfq_taylor([a, b], [c], z0, taylorN)

    F(zn, fn) = ([fn[2], a * b * fn[1] - (c - (a + b + 1) * zn) * fn[2]],
                 [1, zn * (1 - zn)])

    f = ODEpathsolve(z0, [y0, dy0], F, z, h, order = order)[1]

    return f
end

function fastf21(a, b, c, z::Number)
    if abs(z) <= 1
        f = F21(a, b, c, z)
    elseif abs(z) >= 1 && abs(z - 1) >= 1
        z = 1 / z
        f = F21exterior()
    else
        z = 1 - 1 / z
        f = F21interior()
    end

    return f
end

function recursive_2f1_taylor(a, b, c, z0::T, c0, N) where {T <: Number}
    coeffs = zeros(T, N + 1)                                # Initialize taylor coefficients
    coeffs[1:2] .= c0                                       # First 2 coefficients

    A(n) = -(n + a) * (n + b)                               # c_n recursive coefficient
    B(n) =  (n + 1) * (n + c - (1 + a + b + 2n) * z0)       # c_(n+1) recursive coefficient
    C(n) =  (n + 1) * (n + 2) * (1 - z0) * z0               # c_(n+2) recursive coefficient

    # Recursively compute N Taylor coefficients for 2F1
    for k ∈ 0:(N - 2)
        coeffs[k + 3] = -(A(k) * coeffs[k + 1] + B(k) * coeffs[k + 2]) / C(k)
    end

    return Polynomial(coeffs)
end

function fast_2f1(a, b, c, z; z0 = 0im, H = 0.1, N = 150, order = 20)
    if iszero(z0)
        z0 = sgn(imag(z)) * 0.5im
    end

    dir = sign(z - z0)

    zn = z0
    fn = pfq_taylor([a, b], [c], z0, N)

    for i ∈ 1 : ceil(Int64, abs(z - z0) / H)
        h = dir * min(H, abs(z - zn))

        Tn = recursive_2f1_taylor(a, b, c, zn, fn, order)
        Tnp = derivative(Tn)

        zn += h
        fn = [Tn(h), Tnp(h)]
    end

    return fn[1]
end
