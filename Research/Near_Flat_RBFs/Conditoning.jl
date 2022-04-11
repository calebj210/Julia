#=
# Experimenting with the condition number of collocation matrix
#
# Author: Caleb Jacobs
# DLM: 08-04-2022
=#

using Plots
using LinearAlgebra

Φ(r, ε) = ℯ^(-(ε * r)^2)

"""
    colloc(X, ϕ, ε)

Construct the collocation matrix given nodes `X`, radial function `ϕ`, and shape paremeter `ε`.
"""
function colloc(X, ϕ, ε)
    n = size(X, 1)

    return [ϕ(norm(X[:,i] - X[:,j], ε)) for i ∈ 1:n, j ∈ 1:n]
end

"""
    getConds(X, ϕ, ε)

Compute the condition numbers of a the collocation matrix at various values of the shape paremeter.
"""
function getConds(X, ϕ, ε)
    n = size(data, 1)

    conds = zeros(Float64, n)

    for i ∈ 1:n
        conds[i] = cond(colloc(X, ϕ, ε[i]))
    end

    return conds
end


function driver(n = 25, κ = 1)
    ε = rand(ComplexF64, 100)

    x = -1.0 : 2.0/(n + 1) : 1.0
    y = sqrt.(κ^(-2) .- x.^2) .- κ^(-1)

    X = [x'; y']

    scatter(X[1,:], X[2,:])
end
