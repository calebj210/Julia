#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# DLM: May 17, 2024
=#

include("Grid.jl")
using SpecialFunctions

"""
    getModifiedVand(a, b, z)
Compute the modified Vanermonde matrix for computing pFq expansion weights about z = 1.
"""
function getModifiedVand(a, b, z)
    Aα = [(z - 1)^(i - 1) 
        for z ∈ z, i ∈ 1 : length(z) / 2] # Regular shifted Vandermonde for α coefficients

    Aβ = [oneMinusZα(z, sum(b) - sum(a)) * (z - 1)^(i - 1) 
        for z ∈ z, i ∈ 1 : length(z) / 2] # Modified Vandermonde for β coefficients

    return [Aα Aβ]
end

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
             for l ∈ 0 : 35, k ∈ 0 : N - 1])
    s = convert(ComplexF64, s)

    return 2im * (-1 + 0im)^(d - c + α) * (1 - z)^(d - c + α) * sinpi(α) * gamma(d) / z^(d - 1) * s 

#     γ = a[2]
#     return convert(ComplexF64, -2im * sin(π * γ) * (-1 + 0im)^(-c + d - γ) * (1 - z)^(-c + d - γ) * gamma(d) / z^(d - 1) * 
#             sum([(z - 1)^l * gamma(1 + l - γ) / (factorial(big(l)) * gamma(c - l) * gamma(1 + d - c + l - γ)) for l ∈ 0 : 30]))

end
Φ(a, b, ωa, z::Vector) = [Φ(a, b, ωa, z) for z ∈ z]

"""
    getZ1ExpansionWeights(a, b, z)
Compute z = 1 pFq expansion weights.
"""
# function getZ1ExpansionWeights(a, b, z, f)
#     A = BigFloat.(getModifiedVand(a, b, z))
#     A = getModifiedVand(a, b, z)
# 
#     display(cond(A))
# 
#     ω = A \ BigFloat.(f)
#     ω = A \ f
# 
#     ωα = ω[1 : round(Int64, length(ω) / 2)]
#     ωβ = ω[round(Int64, length(ω) / 2) + 1 : end]
# 
#     return (ωα, ωβ)
# end
function getZ1ExpansionWeights(a, b, ωa, z, f)
#     A = lu([(1 - z)^k for z ∈ z, k ∈ 0 : length(z) - 1])
    A = [(1 - z)^k for z ∈ z, k ∈ 0 : length(z) - 1]

    α = sum(b) - sum(a)                                                 # Branch exponent

    ϕ = Φ(a, b, ωa, z)                                                  # Branch correction

    ba = ϕ ./ (1 .- z).^α / (cispi(2α) - 1)                             # Singular RHS
    bb = f - ϕ / (cispi(2α) - 1)                                        # Regular RHS

    ωa = A \ ba                                                         # Singular weights
    ωb = A \ bb                                                         # Regular weights
    
    display(ba)
    display(ωa)

    display(bb)
    display(ωb)

    return(ωa, ωb)
end

"""
    z1PFQ(a, b, ωa, ωb, z)
"""
function z1PFQ(a::Vector, b::Vector, ωa::Vector, ωb::Vector, z::Number)
#     γ = oneMinusZα(z, sum(b) - sum(a))

    α = sum(b) - sum(a)

#     f = sum(((1 - z)^α * ωa + ωb) .* (1 - z).^(0 : length(ωa) - 1))
    N = length(ωa)
    f = (1 - z)^α * sum(ωa .* ((1 - z).^(0 : N - 1))) #+ sum(ωb .* ((1 - z).^(0 : N - 1)))

#     reg = ((gamma(α) * gamma(1 + α) * gamma(b[1])) / (gamma(b[1] - a[1]) * gamma(b[1] - a[2]))) * 
#     (sum([(gamma(a[1] + j) * gamma(a[2] + j) * (1 - z)^j) / (gamma(a[1]) * gamma(a[2]) * factorial(big(j)) * gamma(j + 1 - α)) for j ∈ 0 : 30]))

#     f += reg

    return f
end
z1PFQ(a, b, ωa, ωb, z::Vector) = [z1PFQ(a, b, ωa, ωb, z) for z ∈ z]
