#=
# ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: June 10, 2025
=#

using MathLink
using ArbNumerics: gamma, hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat
import SpecialFunctions.gamma

include("Initialization.jl")
include("Taylor.jl")

function sgn(x)
    iszero(x) ? -one(x) : sign(x)
end

mathematica_2f1(a, b, c, z) = 
    try 
        Complex(weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z).args...)
    catch 
        Real(weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z))
    end
    
mathematica_pfq(a, b, z) = 
    try 
        Complex(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z).args...)
    catch 
        Real(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z))
    end

johansson_2f1(a, b, c, z; bits = 512)::ComplexF64 = arb_2f1(ArbComplex.((a, b, c, z), bits = bits)...)


function taylor_2f1(a, b, c, z::Number; N = 1000, order = 1000, step_max = Inf, init_max = .5)
    # init_max = min(3abs(c / (a * b)), init_max)
    # if abs(z) <= init_max
    #     return maclaurin_2f1(a, b, c, z, N)[1]
    # end
    #
    # z0 = init_max * sign(z)
    #
    # fn = maclaurin_2f1(a, b, c, z0, N)

    # z0, fn = get_initialization(a, b, c, z)

    z0, fn = initialize(a, b, c, z, init_max)

    if z0 == z
        return fn[1]
    end

    # branch = true
    for _ ∈ 1:N
        h_0 = abs(z0) * exp(-2)
        h_1 = abs(z0 - 1) * exp(-2)
        h_f = abs(z0 - z)

        # dir, branch = get_direction(z0, z, fn..., branch)
        dir = sign(z - z0)
        step_size = min(h_1, h_f, h_0, step_max)

        if step_size != h_f
            h = dir * step_size
            fn = recursive_2f1(a, b, c, z0, fn, h, order)
            z0 += h
        else
            h = z - z0
            return recursive_2f1(a, b, c, z0, fn, h, order)[1]
        end
    end

    return fn[1]
end

# function _2f1(a, b, c, z::Number; step_max = Inf, N = 1000, order = 1000)
#     if real(z) >= 0.5 && abs(1 - z) <= 1
#         if isinteger(c - a - b)
#             f = int_abc_2f1(a, b, c, z)
#             f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
#         else
#             f1 = taylor_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z, step_max = step_max, N = N, order = order)
#             f2 = taylor_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z, step_max = step_max, N = N, order = order)
#
#             g1 = gamma(c) / gamma(c - a) * gamma(c - a - b) / gamma(c - b) * z^(-a)
#             g2 = gamma(c) / gamma(a)     * gamma(a + b - c) / gamma(b)     * (1 - z)^(c - a - b) * z^(a - c)
#
#             f = g1 * f1 + g2 * f2
#         end
#     elseif abs(z) >= 1 && abs(1 - z) >= 1
#     # if abs(z) >= 1 && abs(1 - z) >= 1
#         if isinteger(b - a)
#             f = int_ab_2f1(a, b, c, z)
#             f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
#         else
#             f1 = taylor_2f1(a, a - c + 1,  a - b + 1, 1 / z, step_max = step_max, N = N, order = order)
#             f2 = taylor_2f1(b, b - c + 1, -a + b + 1, 1 / z, step_max = step_max, N = N, order = order)
#
#             g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
#             g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)
#
#             f = g1 * f1 + g2 * f2
#         end
#     else
#         f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
#     end
#
#     return f
# end
