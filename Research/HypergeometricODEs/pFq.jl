#=
#   ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: September 10, 2024
=#

using Polynomials
using SpecialFunctions, ArbNumerics
import SpecialFunctions.gamma
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

gamma(z::Complex{BigFloat}) = gamma(ArbComplex(z))
function fast2f1(a, b, c, z::Number; H = 0.1, N = 150, order = 20)
    a,b,c,z = big.((a,b,c,z))

    if real(z) <= 0.5 && abs(z) <= 1
        f = fast_2f1(a, b, c, z, H = H, N = N, order = order)
    elseif abs(z) >= 1 && abs(z - 1) >= 1
        f1 = fast_2f1(a, a - c + 1,  a - b + 1, 1 / z, H = H, N = N, order = order)
        f2 = fast_2f1(b, b - c + 1, -a + b + 1, 1 / z, H = H, N = N, order = order)

        g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
        g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)

        f = g1 * f1 + g2 * f2
    else
        f1 = fast_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z, H = H, N = N, order = order)
        f2 = fast_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z, H = H, N = N, order = order)

        g1 = gamma(c) / gamma(c - a) * gamma(c - a - b) / gamma(c - b) * z^(-a)
        g2 = gamma(c) / gamma(a) * gamma(a + b - c) / gamma(b) * (1 - z)^(c - a - b) * z^(a - c)

        f = g1 * f1 + g2 * f2
    end

    return f
end

function recursive_2f1_taylor(a, b, c, z0::T, c0, N) where {T <: Number}
    coeffs = zeros(T, N + 1)                                # Initialize Taylor coefficients
    coeffs[1:2] .= c0                                       # First 2 coefficients

    # 10 flop optimization for 3-term recurrence
    a0 = -a * b; a1 =  1 - a - b
    b0 = b1 = c - (1 + a + b) * z0; b2 = 2 - 4z0
    c0 = c1 = c2 = 2z0 * (z0 - 1)

    for n = 3 : N + 1
        # Compute next coefficient
        coeffs[n] = (a0 * coeffs[n - 2] + b0 * coeffs[n - 1]) / c0

        # Update recurrence values
        a1 -= 2;  a0 += a1
        b1 += b2; b0 += b1
        c1 += c2; c0 += c1
    end

    return Polynomial(coeffs)
end

function fast_2f1(a, b, c, z; H = 0.1, N = 150, order = 20)
    if abs(z) <= .3
        return pfq_taylor([a, b], [c], z, N)[1]
    end

    z0 = sign(z) * 0.3im
    dir = sign(z - z0)

    zn = z0
    fn = pfq_taylor([a, b], [c], z0, N)

    for i ∈ 1 : ceil(Int64, abs(z - z0) / H)
        h = dir * min(H, abs(z - zn))

        Tn = recursive_2f1_taylor(a, b, c, zn, fn, order + 1)   # Order increase so derivative hits the desired order
        Tnp = derivative(Tn)

        zn += h
        fn = [Tn(h), Tnp(h)]
    end

    return fn[1]
end
