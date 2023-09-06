#=
# Routines for computing hypergeometric functions
#
# Author: Caleb Jacobs
# DLM: September 5, 2023
=#

using SpecialFunctions
include("Gregory.jl")

"""
    _1F_1(a, b, z)
Compute the ₁F₁(a; b; z) hypergeometric function using rescaled Euler formula.
"""
function _1F_1(a, b, z; n = 40, N = 20)
    α = 1 - a
    β = 1 + a - b

    val = dftInt(exp, 0, z, α = α, β = β, n = n, N = N)

    return val / gamma(a) / gamma(b - a) * gamma(b) * z^(1 - b)
end
"""
    _2F_1(a, b, c, z)
Compute the ₂F₁(a, b; c; z) hypergeometric function using rescaled Euler formula.
"""
function _2F_1(a, b, c, z; n = 40, N = 20)
    α = 1 - b
    β = 1 + b - c

    f(z) = (1 - z)^(-a)

    val = dftInt(f, 0, z, α = α, β = β, n = n, N = N)

    return val / gamma(b) / gamma(c - b) * gamma(c) * z^(1 - c)
end

"""
    taylorpFq(A, B, z; N = 50)
Compute pFq using MacLaurin series of pFq using `N` terms in the series.
"""
function taylorpFq(A, B, z; N = 50)
    sn = an = 1           # Initialize partial sum and sequence term
    a = copy(A)
    b = copy(B)

    for n = 1 : N - 1
        an *= prod(a) / prod(b) / n * z

        a .+= 1
        b .+= 1

        sn += an
    end
        
    return sn
end

"""
    pFq(a⃗, b⃗, z)
Compute the generalized pFq hypergeometric function.
"""
function pFq(a, b, z; n = 30, N = 20, start = true, R = .5, TN = 50)
    if abs(z) <= R
        return taylorpFq(a,b,z, N = TN)         # Use Taylor series if z is too close to the origin
    end

    α = 1 - a[end]                              # Base point singularity order
    β = 1 + a[end] - b[end]                     # Evaluation point singularity order

    # Check for known base cases
    if length(a) == 1 && length(b) == 1
        val = dftInt(exp, 0, z, α = α, β = β, n = n, N = N)                         # 1F1
    elseif length(a) == 2 && length(b) == 1
        val = dftInt(x -> (1 - x)^(-a[1]), 0, z, α = α, β = β, n = n, N = N)        # 2F1
    else
        val = dftInt(x -> pFq(a[1 : end - 1], b[1 : end - 1], x, start = false),    # (p-1)F(q-1)
                     0, z, α = α, β = β, n = n, N = N)
    end

    Γ = z^(1 - b[end]) * gamma(b[end]) / gamma(a[end]) / gamma(b[end] - a[end])     # Constant

    return Γ * val
end

### Okay formula except it uses the rescaled formula for 2F1 and 1F1 when outside of taylor radius
# function pFq(a, b, z; n = 40, N = 20, R = .5)
#     if abs(z) <= R
#         return taylorpFq(a,b,z)
#     end
# 
#     if length(a) == 1 && length(b) == 1
#         return _1F_1(a[1], b[1], z, n = n, N = N)
#     end
# 
#     if length(a) == 2 && length(b) == 1
#         return _2F_1(a[1], a[2], b[1], z, n = n, N = N)
#     end
# 
#     α = 1 - a[end]
#     β = 1 + a[end] - b[end]
#     f(t) = pFq(a[1 : end - 1], b[1 : end - 1], t * z)
#     
#     val = dftInt(f, 0, 1, α = α, β = β, n = n, N = N)
# 
#     Γ = gamma(b[end]) /gamma(a[end]) / gamma(b[end] - a[end])
# 
#     return Γ * val
# end
# Wiki formula
function _pFq(a, b, z; n = 40, N = 20)
    α = 1 - a[end]
    β = 1 + a[end] - b[end]

    if length(a) == 1 && length(b) == 1
        val = dftInt(t -> exp(t * z), 0, 1, α = α, β = β, n = n, N = N)
    elseif length(a) == 2 && length(b) == 1
        val = dftInt(t -> (1 - t * z)^(-a[1]), 0, 1, α = α, β = β, n = n, N = N)
    else
        val = dftInt(t -> pFq(a[1 : end - 1], b[1 : end - 1], t * z),
                     0, 1, α = α, β = β, n = n, N = N)
    end

    Γ = gamma(b[end]) /gamma(a[end]) / gamma(b[end] - a[end])

    return Γ * val
end
