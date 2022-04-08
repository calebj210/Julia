#=
# Experimenting with the condition number of collocation matrix
#
# Author: Caleb Jacobs
# DLM: 08-04-2022
=#

using Plots
using LinearAlgebra

"""
    colloc(X, ϕ, ε)

Construct the collocation matrix given nodes `X`, radial function `ϕ`, and shape paremeter `ε`.
"""
function colloc(X, ϕ, ε)
    n = size(X, 1)

    return [ϕ(norm(X[i] - X[j], ε)) for i ∈ 1:n, j ∈ 1:n]
end
