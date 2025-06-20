#=
# Power iteration method for finding eigenvalues
#
# Author: Caleb Jacobs
# Date last modified: 01-02-2022
=#

using LinearAlgebra
using DoubleFloats

function makeHilbert(n, quad = false)
    # Check if we want quad-precision matrix
    if !quad
        A = [Float64(1) / (i + j - 1) for i = 1 : n, j = 1 : n]
    else
        A = [Double64(1) / (i + j - 1) for i = 1 : n, j = 1 : n]
    end

    return A
end

function powIt(A, maxIts)
    q = z = rand(size(A, 1))    # Initialize eigenvectors
    λ = 0                       # Initialize eigenvalue
    for k = 1 : maxIts
        z = A * q               # Compute next eigenvector guess
        q = z / norm(z)         # Normalize eigenvector
        λ = q' * A * q          # Compute next eigenvalue guess
    end

    return (λ, q)
end

function invIt(A, maxIts)
    q = z = ones(Float64, size(A, 1))     # Initialize eigenvectors
    λ = 0                                 # Initialize eigenvalue

    # Perform inverse iteration to find smallest eigenvalue
    for k = 1 : maxIts
        z = A \ q               # Compute next eigenvector guess
        q = z / norm(z)         # Normalize eigenvector
        λ = dot(q, A, q)        # Compute next eigenvalue guess
    end

    return λ
end

function prob2()
    a = 1;          # Lower size of matrix
    b = 10;         # Upper size of matrix
    maxIts = 100;   # Number of iterations for powIt

    data = Vector{Tuple{Float64, Vector{Float64}}}(undef, b - a + 1)

    n = [a : b...]  #

    for i ∈ 1 : b - a + 1
        A = makeHilbert(n[i])

        data[i] = powIt(A, maxIts)
    end

    for i ∈ data
        print(round.(i[1]; digits = 2))
        println(round.(i[2]; digits = 2))
    end
    return(data)
end

function prob3()
    AD = makeHilbert(16)        # Double precision hilbert
    AQ = makeHilbert(16, true)  # Quad precision hilbert

    doubVal = invIt(AD, 1000)   # Double precision λ
    quadVal = invIt(AQ, 1000)   # Quad precision λ

    return (doubVal, quadVal)
end
