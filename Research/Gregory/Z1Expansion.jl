#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# DLM: June 20, 2024
=#

include("Grid.jl")
using SpecialFunctions
using FFTW

binomΓ(n, k) = gamma(n + 1) / gamma(k + 1) / gamma(n - k + 1)

"""
    getZ1ExpansionWeights(a, b, z)
Compute z = 1 pFq expansion weights.
"""
function getZ1ExpansionWeights(a, b, ωa, z, f; n = 150)
    c = a[1]                                                    # Current a coefficient
    d = b[1]                                                    # Current b coefficient
    γ = sum(b) - sum(a) + c - d                                 # Branch exponent

    
    Ak = [π * csc(π * (c - γ - d)) * (d - c) * binomΓ(d - 1, c - 1) *
            sum([(-1)^j * ωa[j + 1] * binomΓ(d + γ + k - 1, d + γ + j - 1) * 
                                      binomΓ(d - c + k - j - 1, -γ - j - 1) 
            for j ∈ 0 : min(k, length(ωa) - 1)]) for k ∈ 0 : n - 1]

    sing = [oneMinusZα(zi, γ - c + d) * sum(Ak .* (1 - zi).^(0 : length(Ak) - 1)) for zi ∈ z]

    A = [(1 - zi)^k for zi ∈ z, k ∈ 0 : n - 1]
    Bk = A \ (f - sing)
    
    return(Ak, Bk)
end

"""
    z1PFQ(a, b, ωa, ωb, z; branch = false)
"""
function z1PFQ(a::Vector, b::Vector, ωa::Vector, ωb::Vector, z::Number; branch = false)
    α = sum(b) - sum(a)
    γ = oneMinusZα(z, α, branch)
    zi = (1 - z).^(0 : length(ωa) - 1)

    f = sum((γ * ωa + ωb) .* zi)
    f = sum(γ * ωa .* zi) + sum(ωb .* zi)

    return f
end
z1PFQ(a, b, ωa, ωb, z::Vector; branch = false) = [z1PFQ(a, b, ωa, ωb, z, branch = branch) for z ∈ z]
