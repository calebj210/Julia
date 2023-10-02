#=
# Testing suite for generalized gregory based quadratures
#
# Author: Caleb Jacobs
# DLM: October 2, 2023
=#

include("GenGreg.jl")

function setupDiff(n, r; α = 0, β = 0, ir = 0.25, np = 3)
    D  = getDiffMat(n, r, α = α, β = β, ir = ir, np = np)   # Differentiation matrix
    g1 = getGrid(n, r, ir = ir, np = p)                     # Padded grid
    g2 = getGrid(n, r, ir = ir, np = 0)                     # Non-padded grid

    return (D, g1, g2)
end

