#=
#   ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: August 30, 2024
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

function fast2f1(a, b, c, z::Number; H = 0.1, N = 150, order = 20)
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

#     A(n) = -(n + a) * (n + b)                               # c_n recursive coefficient
#     B(n) =  (n + 1) * (n + c - (1 + a + b + 2n) * z0)       # c_(n+1) recursive coefficient
#     C(n) =  (n + 1) * (n + 2) * (1 - z0) * z0               # c_(n+2) recursive coefficient
# 
#     # Recursively compute N Taylor coefficients for 2F1
#     for k ∈ 0:(N - 2)
#         coeffs[k + 3] = -(A(k) * coeffs[k + 1] + B(k) * coeffs[k + 2]) / C(k)
#     end

    # 10 flop optimization
    d1 = -a * b
    d2 = 1 - a - b
    d3 = d4 = c - (1 + a + b) * z0;
    d5 = 2 - 4z0
    d6 = d7 = d8 = 2z0 * (z0 - 1)

    for n = 3 : N + 1
        coeffs[n] = (d1 * coeffs[n - 2] + d3 * coeffs[n - 1]) / d6

        d2 -= 2
        d1 += d2
        d4 += d5
        d3 += d4
        d7 += d8
        d6 += d7
    end

    return Polynomial(coeffs)
end

function fast_2f1(a, b, c, z; H = 0.1, N = 150, order = 20)
    if abs(z) <= .5
        return pfq_taylor([a, b], [c], z, N)[1]
    end

    z0 = sign(z) * 0.5im

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
