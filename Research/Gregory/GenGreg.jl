#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: September 8, 2023
=#

using SpecialFunctions
include("Grid.jl")

" Generate vandermond type matrix from nodes in the grid `g` with indices `idx`."
function getVand(idx, g::Grid)
    A = [zk^i for i ∈ 0 : length(idx) - 1, zk ∈ g.z[idx]]

    return A
end

"""
    getIntWeights(idx, g, α, β)

Compute internal weights for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` where `idx` are the indices of the internal nodes in the grid `g`.
"""
function getIntWeights(z, idx, g::Grid, α, β)
    N = length(idx)

    # Left hand side
    A = getVand(idx, g)

    # Right hand side from Taylor expansion of integral
    b = [z^(1 + α + β + k) * gamma(1 + α + k) / gamma(2 + α + β) * gamma(1 + β) for k ∈ 0 : N - 1]

    ω = A \ b   
    
    return ω
end 

"""
    getEndCorrectionWeights(idx, g, α)

Compute end correction weights for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` with singularity of order `α` with points on the grid `g` with indices `idx`.
"""
function getEndCOrrectionWeights(idx, g::Grid, α)
    N = length(idx)                             # Number of correction nodesj

    A = getVand(idx, g)                         # Left hand side

    b = [-zeta(-α - k) / k! for k ∈ 0 : N - 1]  # Right hand side

    ω = A \ b                                   # Solve for weights

    ω *= g.h^(1 + α)                            # Scale weights by grid spacing factor

    return ω
end

"""
    getDiffMat(n, r, α, β)

Generate differentiation matrix for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` over a grid of radius `r`.
"""
function getDiffMat(n, r, α, β)
    
end
