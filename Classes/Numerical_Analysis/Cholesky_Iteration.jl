#=
# Cholesky iteration for finding eigenvalues
#
# Author: Caleb Jacobs
# Date last modified: 20-02-2022
=#

using LinearAlgebra

function CI(A, maxIts)
    for i ∈ 1 : maxIts
        C = cholesky(A)
        A = C.L' * C.L
    end

    return diag(A)
end

function main()
    L = [1.0 0; 2 3]
    A = L * L'

    CI(A, 100) - sort(eigvals(A), rev = true)
end
