#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# DLM: May 30, 2024
=#

include("Grid.jl")
using SpecialFunctions
using FFTW

function getZVals(;r = .25, n = 20)
    z = 1 .+ r * cispi.(-2range(0, 1 - 1 / n, length = n))
end

function getInterpolantNodes(g::Grid, z; n = 1)
    Nx = round(Int64, real(z) / g.h)
    Ny = round(Int64, imag(z) / g.h)

    zIdx = g.c + Nx * g.dx + Ny * g.dy

#     nIdx = vec([zIdx + i * g.dx + j * g.dy for i ∈ -n : n, j ∈ -n : n])
    nIdx = vec([zIdx + i * g.dy for i ∈ -n : n])

    return nIdx
end

"barycentricWeights(z)"
function barycentricWeights(z::Vector)
    n = length(z)

    ω = [1 / prod(z[j] .- z[1 : n .!= j]) for j ∈ 1 : n]

    return ω
# 
#     ω = zeros(ComplexF64, n, n)
# 
#     ω[1, 1] = 1
#     for j ∈ 2 : n
#         for k ∈ 1 : j - 1
#             ω[k, j] = (z[k] - z[j]) * ω[k, j - 1]
#         end
#         ω[j, j] = prod([z[j] - z[k] for k ∈ 1 : j - 1])
#     end
#     ω[:, n] = 1 ./ ω[:, n]
# 
#     return ω[:, n]
end

"""
    barycentricInterpolate(zⱼ, ω, f)
"""
function barycentricInterpolate(zⱼ::Vector, ω::Vector, f::Vector)
    tmp(z) = ω ./ (z .- zⱼ)

    p(z) = any(z .≈ zⱼ) ? only(f[z .≈ zⱼ]) : sum(tmp(z) .* f) / sum(tmp(z))

    return p
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
             for l ∈ 0 : 170, k ∈ 0 : N - 1])
    s = convert(ComplexF64, s)

    return 2im * (-1 + 0im)^(d - c + α) * oneMinusZα(z, d - c + α) * sinpi(α) * gamma(d) / z^(d - 1) * s 
end
Φ(a, b, ωa, z::Vector) = [Φ(a, b, ωa, z) for z ∈ z]

"""
    getZ1ExpansionWeights(a, b, z)
Compute z = 1 pFq expansion weights.
"""
function getZ1ExpansionWeights(a, b, ωa, z, f)
#     r = abs(1 - z[1])
#     D = diagm((-r).^(0 : length(z) - 1))
#     A = [(1 - z)^k for z ∈ z, k ∈ 0 : length(z) - 1] / D
#     A = [(1 - z)^k for z ∈ z, k ∈ 0 : length(z) - 1]
#     display(A)
#     display(cond(A))

    α = sum(b) - sum(a)                                                 # Branch exponent

    ϕ = Φ(a, b, ωa, z)                                                  # Branch correction

    oMz = [oneMinusZα(z, α) for z ∈ z]
    ba = ϕ ./ oMz / (cispi(2α) - 1)                                     # Singular RHS
    bb = f - ϕ / (cispi(2α) - 1)                                        # Regular RHS

#     ωa =  ((-r).^(-(0 : length(z) - 1))).\ (A \ ba)                                                         # Singular weights
#     ωb = ((-r).^(-(0 : length(z) - 1))) .\ (A \ bb)                                                         # Regular weights
#     ωa = A \ ba                                                         # Singular weights
#     ωb = A \ bb                                                         # Regular weights
    ωa = (-abs(1 - z[1])).^(-(0 : length(z) - 1)) .* ifft(ba)
    ωb = (-abs(1 - z[1])).^(-(0 : length(z) - 1)) .* ifft(bb)

#     display(ωa - ωa1)
#     display(ωb - ωb1)
    
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