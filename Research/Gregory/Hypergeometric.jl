#=
# Routines for computing hypergeometric functions
#
# Author: Caleb Jacobs
# DLM: September 1, 2023
=#

using SpecialFunctions
include("Gregory.jl")

"""
    _1F_1(a, b, z)
Compute the ₁F₁(a; b; z) hypergeometric function.
"""
function _1F_1(a, b, z; n = 40, N = 20)
    α = 1 - a
    β = 1 + a - b

    val = dftInt(exp, 0, z, α = α, β = β, n = n, N = N)

    return val / gamma(a) / gamma(b - a) * gamma(b) * z^(1 - b)
end
"""
    _2F_1(a, b, c, z)
Compute the ₂F₁(a, b; c; z) hypergeometric function.
"""
function _2F_1(a, b, c, z; n = 40, N = 20)
    α = 1 - b
    β = 1 + b - c

    f(z) = (1 - z)^(-a)

    val = dftInt(f, 0, z, α = α, β = β, n = n, N = N)

    return val / gamma(b) / gamma(c - b) * gamma(c) * z^(1 - c)
end

"""
    pFq(a⃗, b⃗, z)
Compute the generalized pFq hypergeometric function.
"""
function pFq(a, b, z; n = 40, N = 20, start = true)
    if length(a) == 1 && length(b) == 1
        return _1F_1(a[1], b[1], z, n = n, N = N)
    end

    if length(a) == 2 && length(b) == 1
        return _2F_1(a[1], a[2], b[1], z, n = n, N = N)
    end

    α = 1 - a[end]
    β = 1 + a[end] - b[end]
    f(t) = pFq(a[1 : end - 1], b[1 : end - 1], t * z)
    
    val = dftInt(f, 0, 1, α = α, β = β, n = n, N = N)

    Γ = gamma(b[end]) /gamma(a[end]) / gamma(b[end] - a[end])

    return Γ * val
end
# function pFq(a, b, z; n = 40, N = 20, start = true)
#     if length(a) == 1 && length(b) == 1
#         return _1F_1(a[1], b[1], z, n = n, N = N)
#     end
# 
#     if length(a) == 2 && length(b) == 1
#         return _2F_1(a[1], a[2], b[1], z, n = n, N = N)
#     end
# 
# 
#     α = 1 - a[end]
#     β = 1 + a[end] - b[end]
#     f(z) = pFq(a[1 : end - 1], b[1 : end - 1], z)
#     
#     val = dftInt(f, 0, z, α = α, β = β, n = n, N = N)
# 
#     Γ = gamma(b[end]) /gamma(a[end]) / gamma(b[end] - a[end])
# 
#     return Γ * z^(2 - b[end]) * val
# end
# function pFq(a, b, z; n = 40, N = 20, start = true)
#     if length(a) == 1 && length(b) == 1
#         return _1F_1(a[1], b[1], z, n = n, N = N)
#     end
# 
#     if length(a) == 2 && length(b) == 1
#         return _2F_1(a[1], a[2], b[1], z, n = n, N = N)
#     end
# 
#     α = b[end - 1] - a[end]
#     β = 1 + a[end] - b[end]
# 
#     if abs(α) >= 1
#         display("α = ")
#         display(α)
#     end
#     if abs(β) >= 1
#         display("β = ")
#         display(β)
#     end
# 
#     if β == -1
#         val = dftInt(x -> (z-x) * pFq(a[1 : end - 1], b[1 : end - 1], x, start = false),
#                      0, z, α = α, n = n, N = N)
#     else 
#         val = dftInt(z -> pFq(a[1 : end - 1], b[1 : end - 1], z, start = false),
#                      0, z, α = α, β = β, n = n, N = N)
#     end
#     display("val = ")
#     display(val)
# 
#     Γ = gamma(b[end]) / gamma(a[end]) / gamma(b[end] - a[end])
# 
#     if start
#         return Γ * z^(1 - b[end]) * val
#     else
#         return Γ * val
#     end
# end
