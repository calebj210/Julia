#=
# 2D Harmonic Finite Difference Generator
#
# Author: Caleb Jacobs
# Date last modified: March 4, 2024
=#

using LinearAlgebra
using GenericLinearAlgebra
using Plots
using TaylorDiff

"""
    getCenteredStencil(n = 3)

Generate stencil nodes for a centered `n` × `n` stencil.
"""
function getCenteredStencil(n = 3)
    if n % 2 == 0
        error("Even stencil size not supported.")
    end

    nL = floor(Int64, n / 2)                    # Lower bound on stencil
    Z = [[x, y] for y ∈ -nL : nL, x ∈ -nL : nL] # Stencil nodes

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
    tmp = zeros(1, length(Z))               # Impose central weight is 1
    tmp[ceil(Int64, length(Z) / 2)] = 1     # 1 constraint

    A = [tmp; getA(Z, N)]                   # Construct LHS

    b = [1; zeros(size(A, 1) - 1)]          # Construct Laplacian RHS

    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end

"""
    getD1Weights(Z, N = 7)

Compute stencil weights for first x-derivative at stencil nodes `Z` up to order `N` if possible.
"""
function getD1Weights(Z, N = 7)
    cent = ceil(Int64, length(Z) / 2)
    zrs = zeros(Int64((size(Z, 1) - 1) / 2), length(Z))
    for i = 0 : size(zrs, 1) - 1
        zrs[i + 1, cent + i] = 1
    end

    A = [getA(Z, N); zrs]                   # Construct LHS

    b = zeros(size(A, 1))
    b[2] = 1

    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end

"""
    getD2Weights(Z, N = 7)

Compute stencil weights for second x-derivative at stencil nodes `Z` up to order `N` if possible.
"""
function getD2Weights(Z, N = 7)
    cent = ceil(Int64, length(Z) / 2)
    zrs = zeros(1, length(Z))
    zrs[cent] = 1

    A = [getA(Z, N); zrs]                   # Construct LHS

    b = zeros(size(A, 1))
    b[end] = 1
    b[3] = 2 

    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end

"""
    getD3Weights(Z, N = 7)

Compute stencil weights for first x-derivative at stencil nodes `Z` up to order `N` if possible.
"""
function getD3Weights(Z, N = 7)
    cent = ceil(Int64, length(Z) / 2)
    zrs = zeros(Int64((size(Z, 1) - 1) / 2), length(Z))
    for i = 0 : size(zrs, 1) - 1
        zrs[i + 1, cent + i] = 1
    end

    A = [getA(Z, N); zrs]                   # Construct LHS

    b = zeros(size(A, 1))
    b[4] = 6 

    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end

"""
    getD4Weights(Z, N = 7)

Compute stencil weights for second x-derivative at stencil nodes `Z` up to order `N` if possible.
"""
function getD4Weights(Z, N = 7)

    A = getA(Z, N)                          # Construct LHS

    b = zeros(size(A, 1))
    b[5] = 24 

    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end

