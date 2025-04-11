#=
# 2D Harmonic Finite Difference Generator
#
# Author: Caleb Jacobs
# DLM: April 11, 2025
=#

using LinearAlgebra
using GenericLinearAlgebra
using TaylorDiff

"""
    getCenteredStencil(n = 3)

Generate stencil nodes for a centered `n` × `n` stencil.
"""
function get_centered_stencil(n = 3)
    if n % 2 == 0
        error("Even stencil size not supported.")
    end

    nL = floor(Int64, n / 2)                        # Lower bound on stencil
    Z = [[x, y] for y ∈ -nL:nL, x ∈ -nL:nL] / nL    # Stencil nodes

    return Z
end

"""
    getA(n, N = 7)

Generate LHS for computing 2D harmonic FD stencils with `n`×`n` centered nodes up to order `N`.
"""
function getA(n = 3, N = 7)
    z = vec(get_centered_stencil(n))

    z = [BigFloat.(z) for z ∈ z]

    r = [norm(z) for z ∈ z]                                 # Polar radius of stencil nodes
    θ = [atan(z[2],z[1]) for z ∈ z]                         # Polar angle of stencil nodes
    
    # Generate cosine and sine basis terms
    A_cos = [r^i * cos(i * θ) for i ∈ 1 : N, (r, θ) ∈ zip(r, θ)]
    A_sin = [r^i * sin(i * θ) for i ∈ 1 : N, (r, θ) ∈ zip(r, θ)]

    A = [ones(1, length(z)); 
         reshape([A_cos[:] A_sin[:]]', 2N, length(z))
        ]                                                   # Interlace cosine and sine terms

    return A
end

"""
    getΔWeights(Z, N = 7)

Compute stencil weights for 2D laplacian at stencil nodes `Z` up to order 'N' if possible.

May need to set higher precision via `setprecision()' to get reliable results.
"""
function getΔWeights(n, N = 7)
    Z = get_centered_stencil(n)

    A = [zeros(1, length(Z)); getA(n, N)]
    A[1, ceil(Int64, length(Z) / 2)] = 1

    b = zeros(size(A, 1))                  # Construct Laplacian RHS
    b[1] = 1

    ω = A \ b                                       # Compute weights

    return (Z, reshape(ω, isqrt(length(ω)), :))     # Return weights in Fornberg stencil form
end

"""
    getD1Weights(Z, N = 7)

Compute stencil weights for first x-derivative at stencil nodes `Z` up to order `N` if possible.
"""
function getD1Weights(n, N = 7)
    Z = get_centered_stencil(n)
    A = getA(n, N)                          # Construct LHS

    b = zeros(size(A, 1))
    b[2] = 1

    ω = A \ b                               # Compute weights
    # ω = pinv(A) * b

    return (Z, reshape(ω, isqrt(length(ω)), :))  # Return weights in Fornberg stencil form
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
