#=
# Power iteration method for finding eigenvalues
#
# Author: Caleb Jacobs
# Date last modified: 30-1-2022
=#

using LinearAlgebra

function makeHilbert(n)
    A = [1 / (i + j - 1) for i = 1:n, j = 1:n]

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

    return λ, q
end

function driver(n, maxIts)
    A = makeHilbert(n)

    λ, λ⃗ = powIt(A, maxIts)

    display(λ)

    display(eigmax(A))

    # display(λ⃗)
end
