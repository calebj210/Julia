#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# DLM: May 30, 2024
=#

include("Grid.jl")
using SpecialFunctions
using FFTW

binomΓ(n, k) = gamma(n + 1) / gamma(k + 1) / gamma(n - k + 1)

function getZVals(;r = .25, n = 20)
    z = 1 .- r * cispi.(-2range(0, 1 - 1 / n, length = n))
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
function Φ(a::Vector, b::Vector, ωa::Vector, z::Number; n = 150)
    N = length(ωa)

    c = a[1]                                                    # Current a coefficient
    d = b[1]                                                    # Current b coefficient
    α = sum(b) - sum(a) + c - d                                 # Branch exponent

    s = sum([(ωa[k + 1] * (-1)^l * (1 - z)^(k + l) * gamma(1 + k + l + α)) / 
             (factorial(big(l)) * gamma(c - l) * gamma(1 + d - c + k + l + α))  
             for l ∈ 0 : n, k ∈ 0 : min(19, N - 1)])
    s = convert(ComplexF64, s)

    return 2im * (-1 + 0im)^(d - c + α) * oneMinusZα(z, d - c + α) * sinpi(α) * gamma(d) / z^(d - 1) * s 
end
Φ(a, b, ωa, z::Vector; n = 150) = [Φ(a, b, ωa, z, n = n) for z ∈ z]

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
            for j ∈ 0 : min(k, length(ωa) - 1)]) for k ∈ 0 : length(z) - 1]

    sing = [oneMinusZα(zi, γ - c + d) * sum(Ak .* (1 - zi).^(0 : length(Ak) - 1)) for zi ∈ z]
    
    r  = abs(1 - z[1])                                          # Radius of correction
    scale = r.^(-(0 : length(z) - 1))                           # Coefficient scaling

    Bk = scale .* ifft(f - sing)                                # Regular coefficients
#     display(Bk)

    return(Ak, Bk)
end
# function getZ1ExpansionWeights(a, b, ωa, z, f; n = 150)
#     α = sum(b) - sum(a)                                                 # Branch exponent
# 
#     ϕ = Φ(a, b, ωa, z, n = n)                                           # Branch correction
# 
#     oMz = [oneMinusZα(z, α) for z ∈ z]
#     ba = ϕ ./ oMz / (cispi(2α) - 1)                                     # Singular RHS
#     bb = f - ϕ / (cispi(2α) - 1)                                        # Regular RHS
# 
#     r = -abs(1 - z[1])
#     ωa = r.^(-(0 : length(z) - 1)) .* ifft(ba)
#     ωb = r.^(-(0 : length(z) - 1)) .* ifft(bb)
# 
#     return(ωa, ωb)
# end

"""
    z1PFQ(a, b, ωa, ωb, z)
"""
function z1PFQ(a::Vector, b::Vector, ωa::Vector, ωb::Vector, z::Number)
    α = sum(b) - sum(a)
    γ = oneMinusZα(z, α)
    zi = (1 - z).^(0 : length(ωa) - 1)

    f = sum((γ * ωa + ωb) .* zi)
    f = sum(γ * ωa .* zi) + sum(ωb .* zi)

    return f
end
z1PFQ(a, b, ωa, ωb, z::Vector) = [z1PFQ(a, b, ωa, ωb, z) for z ∈ z]
