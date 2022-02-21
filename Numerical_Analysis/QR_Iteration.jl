#=
# QR iteration for tridiagonal matrices
#
# Author: Caleb Jacobs
# Date last modified: 19/02/2022
=#

using LinearAlgebra

function QRI(A, maxIts; μ = 0)
    for i ∈ 1 : maxIts
        T = qr(A - μ*I)
        A = T.R * T.Q + μ*I
    end

    return A
end

function prob3(n)
    d  = rand(n)
    dl = rand(n - 1)
    A  = Symmetric(Tridiagonal(dl, d, dl))
    B = QRI(A, 200)
    trueλ   = sort(eigvals(Matrix(A)), rev = true)
    approxλ = sort(diag(B), rev = true)

    display(trueλ)
    display(approxλ)

    display(norm(trueλ - approxλ, Inf))
end
