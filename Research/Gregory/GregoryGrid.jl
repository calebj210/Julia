#=
# Generalized Gregory quadrature for computing hypergeometric pFq on a grid
#
# Author: Caleb Jacobs
# DLM: September 6, 2023
=#

"""
    getGrid(n, r)
Generate a complex grid of radius `r` with x and y grid spacing of `n`.
"""
function getGrid(n, r)
    x⃗ = [range(-r, r, length = n)...]
    y⃗ = [range(-r, r, length = n)...]

    grid = x⃗' .+ im * y⃗

    return grid
end
