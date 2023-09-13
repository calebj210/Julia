#=
# Testing suite for generalized gregory based quadratures
#
# Author: Caleb Jacobs
# DLM: September 12, 2023
=#

include("GenGreg.jl")

function setupDiff(n, r; α = 0, β = 0, ir = 0.5, p = 3)
    D = getDiffMat(n, r, α = α, β = β, ir = ir, er = p) # Differentiation matrix
    g1 = getGrid(n, r, ir = ir, p = p)                  # Padded grid
    g2 = getGrid(n, r, ir = ir, p = 0)                  # Non-padded grid

    return (D, g1, g2)
end

