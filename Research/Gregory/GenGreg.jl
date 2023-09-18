#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: September 12, 2023
=#

using SpecialFunctions
include("Grid.jl")

struct Corrections
    bpu::Vector{ComplexF64}
    bpd::Vector{ComplexF64}
    bpl::Vector{ComplexF64}
    bpr::Vector{ComplexF64}

    cru::Vector{ComplexF64}
    crd::Vector{ComplexF64}
    crl::Vector{ComplexF64}
    crr::Vector{ComplexF64}

    epu::Vector{ComplexF64}
    epd::Vector{ComplexF64}
    epl::Vector{ComplexF64}
    epr::Vector{ComplexF64}
end

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
    b = [z^(1 + α + β + k) * gamma(1 + α + k) * gamma(1 + β) / gamma(2 + α + β + k) for k ∈ 0 : N - 1]

    ω = A \ b   
    
    return ω
end 

"""
    getCorrectionWeights(idx, g, α)

Compute end correction weights for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` with singularity of order `α` with points on the grid `g` about the square radius `r`.
"""
function getCorrectionWeights(r, g::Grid; α = 0, dir = 1 + 0im)
    idx = getCorrectionIndices(g.c, g)                      # Correction indices about the origin

    N = length(idx)                                         # Number of correction nodes

    A = getVand(idx, g)                                     # Left hand side

    b = [-zeta(-α - k) * (dir * g.h)^k for k ∈ 0 : N - 1]   # Right hand side

    ω = A \ b                                               # Solve for weights

    ω *=  (dir * g.h)^(1 + α)                               # Scale weights by grid spacing factor

    return ω
end

"""
    getCorrections(g; α = 0, β = 0)

Compute all possible end/corner corrections over a grid `g`.
"""
function getCorrections(g::Grid; α = 0, β = 0)
    # Basepoint corrections
    bpu = getCorrectionWeights(g.np, g, α = α, dir =  0 + 1im)  # Up
    bpd = getCorrectionWeights(g.np, g, α = α, dir =  0 - 1im)  # Down
    bpl = getCorrectionWeights(g.np, g, α = α, dir = -1 + 0im)  # Left
    bpr = getCorrectionWeights(g.np, g, α = α, dir =  1 + 0im)  # Right

    # Corner corrections
    cru = getCorrectionWeights(g.np, g, dir =  0 + 1im)         # Up
    crd = getCorrectionWeights(g.np, g, dir =  0 - 1im)         # Down
    crl = getCorrectionWeights(g.np, g, dir = -1 + 0im)         # Left
    crr = getCorrectionWeights(g.np, g, dir =  1 + 0im)         # Right

    # Endpoint corrections
    epu = getCorrectionWeights(g.np, g, α = β, dir =  0 + 1im)  # Up
    epd = getCorrectionWeights(g.np, g, α = β, dir =  0 - 1im)  # Down
    epl = getCorrectionWeights(g.np, g, α = β, dir = -1 + 0im)  # Left
    epr = getCorrectionWeights(g.np, g, α = β, dir =  1 + 0im)  # Right

    return Corrections(bpu, bpd, bpl, bpr, cru, crd, crl, crr, epu, epd, epl, epr)
end

function getCorrection(c::Corrections, dir, type::String)
    if type == "cr"
        if dir == 0.0 + 1im
            return c.cru
        elseif dir == 0.0 - 1im
            return c.crd
        elseif dir == -1.0 - 0im
            return c.crl
        elseif dir == 1.0 - 0im
            return c.crr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
    elseif type == "bp"
        if dir == 0.0 + 1im
            return c.bpu
        elseif dir == 0.0 - 1im
            return c.bpd
        elseif dir == -1.0 - 0im
            return c.bpl
        elseif dir == 1.0 - 0im
            return c.bpr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
    
    elseif type == "ep"
        if dir == 0.0 + 1im
            return c.epu
        elseif dir == 0.0 - 1im
            return c.epd
        elseif dir == -1.0 - 0im
            return c.epl
        elseif dir == 1.0 - 0im
            return c.epr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
        
    else
        throw(DomainError(type, "needs a valid point type of bp, cr, or ep."))
    end
end

"""
    populateExternalRow!(row, idx, path, g::Grid, α, β)

Populate differentiation entries of `row` corresponding to `g`.z[`zIdx`].
"""
function populateExternalRow!(row, zIdx, g::Grid, α, β)
    path = getPath(zIdx, g, g.np)                                   # Get path from origin to node
    N = length(path)                                                # Number of paths to travel
    c = getCorrections(g, α = α, β = β)                             # End/corner corrections

    for (n, p) ∈ pairs(path)
        dir = sign(p.f - p.i)                                       # Compute direction of travel

        αβt = g.z[p.p].^α .* (g.z[zIdx] .- g.z[p.p]).^β             # Trapezoidal αβ factor

        row[p.p] .= dir * g.h * αβt                                 # Compute trapezoidal weights

        iIdx = getCorrectionIndices(p.i, g)                         # Initial point correction indices
        fIdx = getCorrectionIndices(p.f, g)                         # Final point correction indices

        if n == 1
            βt = (g.z[zIdx] .- g.z[iIdx]).^β                        # Basepoint β factor

            row[iIdx] +=  getCorrection(c, dir, "bp") .* βt         # Basepoint corrections
        else
            αβt = g.z[iIdx].^α .* (g.z[zIdx] .- g.z[iIdx]).^β       # Corner αβ factor

            row[iIdx] +=  getCorrection(c, dir, "cr") .* αβt        # Corner corrections
        end

        if n == N
            αt = g.z[fIdx].^α                                       # Endpoint α factor
            
            row[fIdx] += getCorrection(c, -dir, "ep") .* αt         # Endpoint corrections
        else
            αβt = g.z[fIdx].^α .* (g.z[zIdx] .- g.z[fIdx]).^β       # Corner αβ factor

            row[fIdx] += getCorrection(c, -dir, "cr") .* αβt        # Corner corrections
        end
    end
end

"""
    getDiffMat(n, r, α, β; ir = 0.5, er = 5)

Generate differentiation matrix for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` over a grid of radius `r`.

The radius to use the Taylor expansion is given by `ir` while the relative radius of the correction stencils are given by `er`.
"""
function getDiffMat(n, r; α = 0.0, β = 0.0, ir = 0.5, er = 3)
    g = getGrid(n, r, ir = ir, p = er)                                  # Generate grid

    (iMap, eMap) = getReducedGridMap(g)                                 # Index maps to reduced grid for indexing through diff matrix

    D = zeros(ComplexF64, length(g.i) + length(g.e), length(g.z))       # Initialize differentiation matrix

    # Populate internal weights using Taylor expansion approximation
    for (i, iIdx) ∈ pairs(g.i)
        D[iMap[i], g.i] = getInternalWeights(iIdx, g, α, β)
    end

    # Populate external weights using generalized Gregory quadrature
    for (e, eIdx) ∈ pairs(g.e)
        populateExternalRow!(D[eMap[e], :], eIdx, g, α, β)
    end

    return D
end

#     # Populate external weights using trapezoidal rule and end corrections
#     bpω = getCorrectionWeights(er, g, α)                        # Base point correction
#     crω = getCorrectionWeights(er, g, 0.0)                      # Corner correction
#     epω = getCorrectionWeights(er, g, β)                        # End point correction
# 
#     for (i, eIdx) ∈ pairs(g.e)
#         z = g.z[eIdx]                                           # Current z value
# 
#         # Populate trapezoidal weights
#         (hIdx, vIdx, cIdx) = getPath(eIdx, g)                   # Indices along path of integration
# 
#         αβt = g.z[hIdx].^α .* (z .- g.z[hIdx]).^β               # Horizontal α-β term in integrand
#         if real(z) > 0
#             D[eMap[i], hIdx] .=  g.h * αβt                      # Set trapezoidal weights moving right
#         else
#             D[eMap[i], hIdx] .= -g.h * αβt                      # Set trapezoidal weights moving left
#         end
# 
#         αβt = g.z[vIdx].^α .* (z .- g.z[vIdx]).^β               # Vertical α-β term in integrand
#         if imag(z) > 0
#             D[eMap[i], vIdx] .=  im * g.h * αβt                 # Set trapezoidal weights moving up
#         else
#             D[eMap[i], vIdx] .= -im * g.h * αβt                 # Set trapezoidal weights moving down
#         end
# 
#         # Populate left-right end correction weights
#         bpIdx = getCorrectionIndices(g.c,  er, g)               # Base point correction indices
#         crIdx = getCorrectionIndices(cIdx, er, g)               # Corner correction indices
#         epIdx = getCorrectionIndices(eIdx, er, g)               # Endpoint correction indices
#         
#         if imag(g.z[eIdx]) != 0
#             if imag(z) > 0
#                 rbpIdx = rotCorrection(bpIdx, 1)                # Rotated base point indices
#                 rcrIdx = rotCorrection(crIdx, 3)                # Rotated base point indices
# 
#                 hbp = im^(1 + α)                                # Base point angle correction
#                 hcr = im                                        # Corner  angle correction
#             else 
#                 rbpIdx = rotCorrection(bpIdx, 3)                # Rotated base point indices
#                 rcrIdx = rotCorrection(crIdx, 1)                # Rotated base point indices
# 
#                 hbp = (-im)^(1 + α)                             # Base point angle correction
#                 hcr = (-im)                                     # Corner  angle correction
#             end
# 
#             βt  = (z .- g.z[rbpIdx]).^β                         # Base point β term in integrand
#             αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β       # Corner α-β term in integrand
# 
#             D[eMap[i], rbpIdx] += hbp * bpω .* βt               # Basepoint correction
#             D[eMap[i], rcrIdx] += hcr * crω .* αβt              # Corner left-right correction
#         end
# 
#         if real(g.z[eIdx]) != 0
#             if real(z) > 0
#                 rcrIdx = rotCorrection(crIdx, 0)                # Rotated corner indices
#                 repIdx = rotCorrection(epIdx, 2)                # Rotated endpoint indices
# 
#                 hcr = 1.0                                       # Corner angle correction
#                 hep = 1.0^(1 + β)                               # Endpoint angle correction
#             else
#                 rcrIdx = rotCorrection(crIdx, 2)                # Rotated corner indices
#                 repIdx = rotCorrection(epIdx, 0)                # Rotated endpoint indices
# 
#                 hcr = (-1 + 0im)                                # Corner angle correction
#                 hep = (-1 + 0im)^(1 + β)                        # Endpoint angle correction
#             end
# 
#             αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β       # Corner α-β term in integrand
#             αt  = g.z[repIdx].^α                                # Endpoint α term in integrand
# 
#             D[eMap[i], rcrIdx] += hcr * crω .* αβt              # Corner up-down correction
#             D[eMap[i], repIdx] += hep * epω .* αt               # Endpoint correction
#         end

#         if real(z) != 0
#             if real(z) > 0
#                 rbpIdx = rotCorrection(bpIdx, 0)                # Rotated base point indices
#                 rcrIdx = rotCorrection(crIdx, 2)                # Rotated base point indices
# 
#                 hbp = 1^(1 + α)                                 # Base point angle correction
#                 hcr = 1                                         # Corner  angle correction
#             else 
#                 rbpIdx = rotCorrection(bpIdx, 2)                # Rotated base point indices
#                 rcrIdx = rotCorrection(crIdx, 0)                # Rotated base point indices
# 
#                 hbp = (-1 + 0im)^(1 + α)                        # Base point angle correction
#                 hcr = (-1 + 0im)                                # Corner  angle correction
#             end
# 
#             βt  = (z .- g.z[rbpIdx]).^β                         # Base point β term in integrand
#             αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β       # Corner α-β term in integrand
# 
#             D[eMap[i], rbpIdx] += hbp * bpω .* βt               # Basepoint correction
#             D[eMap[i], rcrIdx] += hcr * crω .* αβt              # Corner left-right correction
#         end
# 
#         if imag(z) != 0
#             if imag(z) > 0
#                 rcrIdx = rotCorrection(crIdx, 1)                # Rotated corner indices
#                 repIdx = rotCorrection(epIdx, 3)                # Rotated endpoint indices
# 
#                 hcr = im                                        # Corner angle correction
#                 hep = im^(1 + β)                                # Endpoint angle correction
#             else
#                 rcrIdx = rotCorrection(crIdx, 3)                # Rotated corner indices
#                 repIdx = rotCorrection(epIdx, 1)                # Rotated endpoint indices
# 
#                 hcr = (-im)                                     # Corner angle correction
#                 hep = (-im)^(1 + β)                             # Endpoint angle correction
#             end
# 
#             αβt = g.z[rcrIdx].^α .* (z .- g.z[rcrIdx]).^β       # Corner α-β term in integrand
#             αt  = g.z[repIdx].^α                                # Endpoint α term in integrand
# 
#             D[eMap[i], rcrIdx] += hcr * crω .* αβt              # Corner up-down correction
#             D[eMap[i], repIdx] += hep * epω .* αt               # Endpoint correction
#         end
#     end
# 
#     return D
# end
