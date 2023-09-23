#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: September 20, 2023
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

"Compute roots given a power α and a branch cut rotation of θ."
function zα(z, α; θ = 0)
    if θ >= 0
        return angle(z) >= θ - π ? z^α : z^α * cispi( 2 * (α % 1))
    else
        return angle(z) <  θ + π ? z^α : z^α * cispi(-2 * (α % 1))
    end
end

"Modified sign function to return 1 when z = 0"
sgn(z) = iszero(z) ? one(z) : sign(z)

"Generate vandermond type matrix from nodes in the grid `g` with indices `idx`."
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

    ω *=  g.h^(1 + α)                                       # Scale weights by grid spacing factor

    return ω
end

"""
    getCorrections(g; α = 0, β = 0)

Compute all possible end/corner corrections over a grid `g`.
"""
function getCorrections(g::Grid; α = 0, β = 0)
    # Basepoint corrections
#     bpu = getCorrectionWeights(g.np, g, α = α, dir =  0 + 1im)  # Up
#     bpd = getCorrectionWeights(g.np, g, α = α, dir =  0 - 1im)  # Down
#     bpl = getCorrectionWeights(g.np, g, α = α, dir = -1 + 0im)  # Left
    bpr = getCorrectionWeights(g.np, g, α = α, dir =  1 + 0im)  # Right
    bpu = rotCorrection(bpr, -1)
    bpl = rotCorrection(bpr, -2)
    bpd = rotCorrection(bpr, -3)

    # Corner corrections
#     cru = getCorrectionWeights(g.np, g, dir =  0 + 1im)         # Up
#     crd = getCorrectionWeights(g.np, g, dir =  0 - 1im)         # Down
#     crl = getCorrectionWeights(g.np, g, dir = -1 + 0im)         # Left
    crr = getCorrectionWeights(g.np, g, dir =  1 + 0im)         # Right
    cru = rotCorrection(crr, -1)
    crl = rotCorrection(crr, -2)
    crd = rotCorrection(crr, -3)

    # Endpoint corrections
#     epu = getCorrectionWeights(g.np, g, α = β, dir =  0 + 1im)  # Up
#     epd = getCorrectionWeights(g.np, g, α = β, dir =  0 - 1im)  # Down
#     epl = getCorrectionWeights(g.np, g, α = β, dir = -1 + 0im)  # Left
    epr = getCorrectionWeights(g.np, g, α = β, dir =  1 + 0im)  # Right
    epu = rotCorrection(epr, -1)
    epl = rotCorrection(epr, -2)
    epd = rotCorrection(epr, -3)

    return Corrections(bpu, bpd, bpl, bpr, cru, crd, crl, crr, epu, epd, epl, epr)
end

function getCorrection(c::Corrections, dir, type::String)
    if type == "cr"
        if     dir ==  0.0 + 1im
            return c.cru
        elseif dir ==  0.0 - 1im
            return c.crd
        elseif dir == -1.0 + 0im
            return c.crl
        elseif dir ==  1.0 + 0im
            return c.crr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
    elseif type == "bp"
        if     dir ==  0.0 + 1im
            return c.bpu
        elseif dir ==  0.0 - 1im
            return c.bpd
        elseif dir == -1.0 + 0im
            return c.bpl
        elseif dir ==  1.0 + 0im
            return c.bpr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
    
    elseif type == "ep"
        if     dir ==  0.0 + 1im
            return c.epu
        elseif dir ==  0.0 - 1im
            return c.epd
        elseif dir == -1.0 + 0im
            return c.epl
        elseif dir ==  1.0 + 0im
            return c.epr
        else
            throw(DomainError(dir, "not a valid direction of grid integration"))
        end
        
    else
        throw(DomainError(type, "needs a valid point type of bp, cr, or ep."))
    end
end

"""
    getExternalWeights(row, zIdx, path, g::Grid, α, β)

Get external differentiation entries corresponding to `g`.z[`zIdx`].
"""
function getExternalWeights(zIdx, g::Grid, α, β)
    path = getPath(zIdx, g, 2g.np)                                      # Get path from origin to node
    N = length(path)                                                    # Number of paths to travel
    c = getCorrections(g, α = α, β = β)                                 # End/corner corrections
    z = g.z[zIdx]                                                       # Current z value
    row = zeros(ComplexF64, length(g.z))                                # Initialize row to populate

    for (n, p) ∈ pairs(path)
        dir = sign(g.z[p.f] - g.z[p.i])                                 # Compute direction of travel

        if real(z) >= 2g.h * g.np
            αβt = g.z[p.p].^α .* (z .- g.z[p.p]).^β                     # Trapezoidal αβ factor
        else
            αβt = zα.(     g.z[p.p], α, θ = sgn(imag(z)) * π) .* 
                  zα.(z .- g.z[p.p], β, θ = sgn(imag(z)) * π/2)        # Trapezoidal αβ factor
        end

        row[p.p] .+= dir * g.h * αβt                                    # Compute trapezoidal weights

        iIdx = getCorrectionIndices(p.i, g)                             # Initial point correction indices
        fIdx = getCorrectionIndices(p.f, g)                             # Final point correction indices

        if n == 1
            βt = zα.(z .- g.z[iIdx], β)                                 # Basepoint β factor
            h = zα(dir, 1 + α, θ = sgn(imag(z)) * π/2)

            row[iIdx] +=  h * βt .* getCorrection(c, dir, "bp")         # Basepoint corrections
        else
            if real(z) <= 2g.h * g.np
                αβt = zα.(     g.z[iIdx], α, θ = sgn(imag(z)) * π) .* 
                      zα.(z .- g.z[iIdx], β, θ = sgn(imag(z)) * π/2)   # Corner αβ factor
            else
                αβt = g.z[iIdx].^α .* (z .- g.z[iIdx]).^β               # Corner αβ factor
            end
            h = dir

            row[iIdx] +=  h * αβt .* getCorrection(c, dir, "cr")        # Corner corrections
        end

        if n == N
            if real(z) <= 2g.h * g.np
                αt = zα.(g.z[fIdx], α, θ = sgn(imag(z)) * π)
            else
                αt = g.z[fIdx].^α                                       # Endpoint α factor
            end
            h = dir^(1 + β)
            
            row[fIdx] += h * αt .* getCorrection(c, -dir, "ep")         # Endpoint corrections
        else
            if real(z) <= 2g.h * g.np
                αβt = zα.(     g.z[fIdx], α, θ = sgn(imag(z)) * π) .* 
                      zα.(z .- g.z[fIdx], β, θ = sgn(imag(z)) * π/2)   # Corner αβ factor
            else
                αβt = g.z[fIdx].^α .* (z .- g.z[fIdx]).^β               # Corner αβ factor
            end
            h = dir

            row[fIdx] += h * αβt .* getCorrection(c, -dir, "cr")        # Corner corrections
        end
    end

    return row
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
        D[eMap[e], :] = getExternalWeights(eIdx, g, α, β)
    end

    return D
end
