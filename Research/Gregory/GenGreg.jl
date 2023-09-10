#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: September 10, 2023
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
function getInternalWeights(zIdx, g::Grid, α, β)
    z = g.z[zIdx]               # Get z value at index 

    N = length(g.i)             # Number of internal nodes

    A = getVand(g.i, g)         # Left hand side

    # Right hand side from Taylor expansion of integral
    b = [z^(1 + α + β + k) * gamma(1 + α + k) / gamma(2 + α + β) * gamma(1 + β) for k ∈ 0 : N - 1]

    ω = A \ b   
    
    return ω
end 

"""
    getCorrectionWeights(idx, g, α)

Compute end correction weights for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` with singularity of order `α` with points on the grid `g` about the square radius `r`.
"""
function getCorrectionWeights(r, g::Grid, α)
    idx = getCorrectionIndices(g.c, r, g)           # Correction indices about the origin

    N = length(idx)                                 # Number of correction nodes

    A = getVand(idx, g)                             # Left hand side

    b = [-zeta(-α - k) * g.h^k for k ∈ 0 : N - 1]   # Right hand side

    ω = A \ b                                       # Solve for weights

    ω *= g.h^(1 + α)                                # Scale weights by grid spacing factor

    return ω
end

"""
    getDiffMat(n, r, α, β; ir = 0.5, er = 5)

Generate differentiation matrix for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` over a grid of radius `r`.

The radius to use the Taylor expansion is given by `ir` while the relative radius of the correction stencils are given by `er`.
"""
function getDiffMat(n, r, α, β; ir = 0.5, er = 3)
    g = getGrid(n, r, ir = ir, p = er)                      # Generate grid

    (iMap, eMap) = getReducedGridMap(g)                     # Index maps to reduced grid for indexing through diff matrix

    D = zeros(ComplexF64, length(g.i) + length(g.e), length(g.z))       # Initialize differentiation matrix

    # Populate internal weights using Taylor expansion approximation
    for (i, iIdx) ∈ pairs(g.i)
        D[iMap[i], g.i] = getInternalWeights(iIdx, g, α, β)
    end

    # Populate external weights using trapezoidal rule and end corrections
    bpω = getCorrectionWeights(er, g, α)                    # Base point correction
    crω = getCorrectionWeights(er, g, 0.0)                  # Corner correction
    epω = getCorrectionWeights(er, g, β)                    # End point correction

    for (i, eIdx) ∈ pairs(g.e)
        z = g.z[eIdx]                                       # Current z value

        # Populate trapezoidal weights
        (hIdx, vIdx) = getPathIndices(eIdx, g)              # Indices along path of integration

        αβt = g.z[hIdx].^α .* (z .- g.z[hIdx]).^β           # Horizontal α-β term in integrand
        if real(z) > 0
            D[eMap[i], hIdx] =  g.h * αβt                   # Set trapezoidal weights moving right
        else
            D[eMap[i], hIdx] = -g.h * αβt                   # Set trapezoidal weights moving left
        end

        αβt = g.z[vIdx].^α .* (z .- g.z[vIdx]).^β           # Vertical α-β term in integrand
        if imag(z) > 0
            D[eMap[i], vIdx] .=  im * g.h * αβt             # Set trapezoidal weights moving up
        else
            D[eMap[i], vIdx] .= -im * g.h * αβt             # Set trapezoidal weights moving down
        end

        # Populate left-right end correction weights
        bpIdx = getCorrectionIndices(eIdx, er, g)           # Base point correction indices
        epIdx = getCorrectionIndices(eIdx, er, g)           # Endpoint correction indices
        crIdx = getCorrectionIndices(eIdx, er, g)           # Corner correction indices
        
        if !isempty(hIdx)
            if real(z) > 0
                rcrIdx = rotCorrection(crIdx, 2, g)             # Rotated base point indices

                βt  = (z .- g.z[bpIdx]).^β                      # Base point β term in integrand
                αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β   # Corner α-β term in integrand

                D[eMap[i], bpIdx]  += bpω .* (z .- g.z[bpIdx]).^β
                D[eMap[i], rcrIdx] -= crω .* αβt
            else 
                rbpIdx = rotCorrection(bpIdx, 2, g)             # Rotated base point indices

                βt  = (z .- g.z[rbpIdx]).^β                     # Base point β term in integrand
                αβt = g.z[crIdx].^α .* (z .- g.z[crIdx]).^β     # Corner α-β term in integrand

                D[eMap[i], rbpIdx]  -= bpω .* βt                # Basepoint correction 
                D[eMap[i], crIdx] += crω .* αβt                 # Corner left-right corner correction
            end
        end

        if !isempty(vIdx)
            if imag(z) > 0
                rcrIdx = rotCorrection(crIdx, 1, g)             # Rotated corner indices
                repIdx = rotCorrection(epIdx, 3, g)             # Rotated endpoint indices

                αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β   # Corner α-β term in integrand
                αt  = (z .- g.z[repIdx]).^α                     # Endpoint α term in integrand
                
                D[eMap[i], rcrIdx] += im * crω .* αβt           # Corner up-down correction
                D[eMap[i], repIdx] -= im * epω .* αt            # Endpoint correction
            else
                rcrIdx = rotCorrection(crIdx, 3, g)             # Rotated corner indices
                repIdx = rotCorrection(epIdx, 1, g)             # Rotated endpoint indices

                αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β   # Corner α-β term in integrand
                αt  = (z .- g.z[repIdx]).^α                     # Endpoint α term  in integrand 
                
                D[eMap[i], rcrIdx] -= im * crω .* αβt           # Corner up-down correction
                D[eMap[i], repIdx] += im * epω .* αt            # Endpoint correction
            end
        end
    end

    return D
end
