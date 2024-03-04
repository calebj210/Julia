#=
# 2D Harmonic Finite Difference Generator
#
# Author: Caleb Jacobs
# Date last modified: March 4, 2024
=#

using LinearAlgebra
using GenericLinearAlgebra
using Plots
using ForwardDiff

"""
    getCenteredStencil(n = 3)

Generate stencil nodes for a centered `n` × `n` stencil.
"""
function getCenteredStencil(n = 3)
    if n % 2 == 0
        error("Even stencil size not supported.")
    end

    nL = floor(Int64, n / 2)                    # Lower bound on stencil
    Z = [[x, y] for x ∈ -nL : nL, y ∈ -nL : nL] # Stencil nodes

    return Z
end

"""
    getA(Z, N = 7)

Generate LHS for computing 2D harmonic FD stencils at nodes `Z` up to order `N`.
"""
function getA(Z, N = 7)
    z = vec(Z)                                  # Vectorized stencil

    z = [BigFloat.(z) for z ∈ z]

    r = [norm(z) for z ∈ z]                     # Polar radius of stencil nodes
    θ = [atan(z[2],z[1]) for z ∈ z]             # Polar angle of stencil nodes
    
    # Generate cosine and sine basis terms
    ACos = [r^i * cos(i * θ) for i ∈ 0 : N, (r, θ) ∈ zip(r, θ)]
    ASin = [r^i * sin(i * θ) for i ∈ 1 : N, (r, θ) ∈ zip(r, θ)]

    A = [ACos; ASin]                            # Full LHS

    return A
end

"""
    getΔWeights(Z, N = 7)

Compute stencil weights for 2D laplacian at stencil nodes `Z` up to order 'N' if possible.

May need to set higher precision via `setprecision()' to get reliable results.
"""
function getΔWeights(Z, N = 7)
    tmp = zeros(1, length(Z))
    tmp[ceil(Int64, length(Z) / 2)] = 1

    A = [tmp; getA(Z, N)]

    b = [1; zeros(size(A, 1) - 1)]

    ω = A \ b

    return reshape(ω, isqrt(length(ω)), :)
end
