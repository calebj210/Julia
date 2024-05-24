#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# DLM: May 23, 2024
=#

include("Grid.jl")
using SpecialFunctions

""" 
    Φ(a, b, ωa, z)
Compute branch correction given previous correction weights `ωa` at `z`.
"""
function Φ(a::Vector, b::Vector, ωa::Vector, z::Number)
    N = length(ωa)

    c = a[1]                                                    # Current a coefficient
    d = b[1]                                                    # Current b coefficient
    α = sum(b) - sum(a) + c - d                                 # Branch exponent

    s = sum([(ωa[k + 1] * (-1)^l * (1 - z)^(k + l) * gamma(1 + k + l + α)) / 
             (factorial(big(l)) * gamma(c - l) * gamma(1 + d - c + k + l + α))  
             for l ∈ 0 : 60, k ∈ 0 : N - 1])
    s = convert(ComplexF64, s)

    return 2im * (-1 + 0im)^(d - c + α) * oneMinusZα(z, d - c + α) * sinpi(α) * gamma(d) / z^(d - 1) * s 
end
Φ(a, b, ωa, z::Vector) = [Φ(a, b, ωa, z) for z ∈ z]

"""
    getZ1ExpansionWeights(a, b, z)
Compute z = 1 pFq expansion weights.
"""
function getZ1ExpansionWeights(a, b, ωa, z, f)
    A = lu([(1 - z)^k for z ∈ z, k ∈ 0 : length(z) - 1])

    α = sum(b) - sum(a)                                                 # Branch exponent

    ϕ = Φ(a, b, ωa, z)                                                  # Branch correction

    oMz = [oneMinusZα(z, α) for z ∈ z]
    ba = ϕ ./ oMz / (cispi(2α) - 1)                                     # Singular RHS
    bb = f - ϕ / (cispi(2α) - 1)                                        # Regular RHS

    ωa = A \ ba                                                         # Singular weights
    ωb = A \ bb                                                         # Regular weights
    
    return(ωa, ωb)
end

"""
    z1PFQ(a, b, ωa, ωb, z)
"""
function z1PFQ(a::Vector, b::Vector, ωa::Vector, ωb::Vector, z::Number)
    α = sum(b) - sum(a)
    γ = oneMinusZα(z, α)

    f = sum((γ * ωa + ωb) .* (1 - z).^(0 : length(ωa) - 1))

    return f
end
z1PFQ(a, b, ωa, ωb, z::Vector) = [z1PFQ(a, b, ωa, ωb, z) for z ∈ z]
