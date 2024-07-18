#=
#   ODE approach to computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: July 3, 2024
=#

using Polynomials
using SpecialFunctions
include("TimeStep.jl")

function pFqTaylor(a, b, z, N)
    coeffs = zeros(typeof(z), N + 1)
    coeffs[1] = 1.0

    for n ∈ 1 : N
        coeffs[n + 1] = coeffs[n] * prod(a .+ (n - 1)) / prod(b .+ (n - 1)) / n 
    end

    p  = Polynomial(coeffs)
    dp = derivative(p)

    return (p(z), dp(z))
end

function F21(a, b, c, z::Number; h = .1, z0 = .25 + 0im, order = 20, TaylorN = 50)
    y0, dy0 = pFqTaylor([a, b], [c], z0, TaylorN)

    F(zn, fn) = ([fn[2], a * b * fn[1] - (c - (a + b + 1) * zn) * fn[2]],
                 [1, zn * (1 - zn)])

    f = odeSolve(z0, [y0, dy0], F, z, h, order = order, store = false)[1]

    return f
end
