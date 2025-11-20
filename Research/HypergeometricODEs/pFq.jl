#=
# ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: November 18, 2025
=#

using MathLink
using MATLAB
using ArbNumerics: gamma, hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat
using HypergeometricFunctions: pFqweniger as weniger_pfq
import SpecialFunctions.gamma

# include("Initialization.jl")
# include("Taylor.jl")

function sgn(x)
    iszero(x) ? -one(x) : sign(x)
end

function mathematica_2f1(a, b, c, z)
    val = weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z)
    try 
        out = Complex(val.args...)
        return out
    catch 
        out = Real(val)
        return out
    end
end

matlab_2f1(a, b, c, z) = mat"hypergeom([$a, $b], [$c], $z)"

mat"addpath('/home/merlin/Documents/Papers/Gauss Hypergeometric/Crespo Code/SecondLink/')"
mat"addpath('/home/merlin/Documents/Papers/Gauss Hypergeometric/Crespo Code/SecondLink/chebfun-master/chebfun-master/')"
mat"warning('off', 'MATLAB:singularMatrix')"
mat"warning('off', 'MATLAB:nearlySingularMatrix')"
function uf_2f1(a, b, c, z)
    val = mat"hypergeom_real($a, $b, $c, 160, $z)"

    if length(val) == 1
        return val[1]
    else
        return NaN + NaN*im
    end
end

johansson_2f1(a, b, c, z; bits = 512)::ComplexF64 = arb_2f1(ArbComplex.((a, b, c, z), bits = bits)...)

"Levin-type factorial"
function weniger_2f1(a, b, c, z::Number) 
    (a,b,c,z) = convert.(ComplexF64, (a,b,c,z))
    return weniger_pfq((a,b), (c,), z)   
end

# function taylor_2f1(a, b, c, z::Number; N = 1000, order = 1000, step_max = Inf, init_max = .5)
#     # Loop path
#     rng = 4
#     if real(z) > -rng
#         z0, fn = initialize(a, b, c, -rng + 0im, init_max)
#     else
#         z0, fn = initialize(a, b, c, z, init_max)
#     end
#
#     # z > 1 branch wall errors
#     # if real(z) > 1 && abs(angle(z)) <= π / 4
#     #     z0, fn = initialize(a, b, c, 1 + sgn(imag(z)) * im, init_max)
#     # else
#     #     z0, fn = initialize(a, b, c, z, init_max)
#     # end
#
#     # Default
#     # z0, fn = initialize(a, b, c, z, init_max)
#
#     if z0 == z
#         return fn[1]
#     end
#
#     for _ ∈ 1:N
#         h_0 = abs(z0) * exp(-2)
#         h_1 = abs(z0 - 1) * exp(-2)
#         h_f = abs(z0 - z)
#
#         dir = get_direction(z0, z)
#         step_size = min(h_1, h_f, h_0, step_max)
#
#         if step_size != h_f
#             h = dir * step_size
#             fn = recursive_2f1(a, b, c, z0, fn, h, order)
#             z0 += h
#         else
#             h = z - z0
#             return recursive_2f1(a, b, c, z0, fn, h, order)[1]
#         end
#     end
#
#     return fn[1]
# end

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
