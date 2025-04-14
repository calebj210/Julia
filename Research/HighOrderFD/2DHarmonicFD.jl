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

function solvable(A, b)
    ker = nullspace(A')
    if isapprox(ker' * b, zeros(size(ker, 2)), atol = sqrt(eps(BigFloat)))
        println("Solvable!")
        # display(ker' * b)
        return true
    else
        println("Not solvable :(")
        # display(ker' * b)
        return false
    end
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

    solvable(A, b)

    ω = A \ b                                       # Compute weights

    return reshape(ω, isqrt(length(ω)), :)     # Return weights in Fornberg stencil form
end

function Dd_stencil(n, N, d)
    A = getA(n, N)
    b = 2N + 1

    b = zeros(2N + 1)
    b[2d] = factorial(d)

    solvable(A, b)
    
    ω = A \ b                               # Compute weights

    return reshape(ω, isqrt(length(ω)), :)  # Return weights in Fornberg stencil form
end
