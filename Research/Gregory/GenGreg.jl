#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: October 3, 2023
=#

using SpecialFunctions
using LinearAlgebra

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
function zα(z, α::Real; θ = 0)
    if θ >= 0
        return angle(z) >= θ - π ? z^α : z^α * cispi( 2 * (α % 1))
    else
        return angle(z) <  θ + π ? z^α : z^α * cispi(-2 * (α % 1))
    end
end

"Rotate branch of complex power"
function zα(z, α::Complex; θ = 0)
    if θ >= 0
        return angle(z) >= θ - π ? z^α : z^α * cispi( 2 * (α.re % 1)) * exp(-2π * α.im)
    else
        return angle(z) <  θ + π ? z^α : z^α * cispi(-2 * (α.re % 1)) * exp( 2π * α.im)
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
function getInternalWeights(zIdx, A, g::Grid, α, β)
    z = g.z[zIdx]               # Get z value at index 

    N = size(A, 1)

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
    idx = getCorrectionIndices(g.c, g.np, g)                # Correction indices about the origin

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
    bpr = getCorrectionWeights(g.np, g, α = α, dir =  1 + 0im)  # Basepoint right
    bpu = rotCorrection(bpr, -1)                                # Basepoint up
    bpl = rotCorrection(bpr, -2)                                # Basepoint left
    bpd = rotCorrection(bpr, -3)                                # Basepoint down

    # Corner corrections
    crr = getCorrectionWeights(g.np, g, dir =  1 + 0im)         # Corner right
    cru = rotCorrection(crr, -1)                                # Corner up
    crl = rotCorrection(crr, -2)                                # Corner left
    crd = rotCorrection(crr, -3)                                # Corner down

    # Endpoint corrections
    epr = getCorrectionWeights(g.np, g, α = β, dir =  1 + 0im)  # Endpoint right
    epu = rotCorrection(epr, -1)                                # Endpoint up
    epl = rotCorrection(epr, -2)                                # Endpoint left
    epd = rotCorrection(epr, -3)                                # Endpoint down

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
function getExternalWeights(zIdx, c::Corrections, g::Grid, α, β)
    path = getPath(zIdx, g, 2g.np)                                      # Get path from origin to node
    N = length(path)                                                    # Number of paths to travel
    z = g.z[zIdx]                                                       # Current z value
    row = zeros(ComplexF64, length(g.z))                                # Initialize row to populate

    for (n, p) ∈ pairs(path)
        dir = sign(g.z[p.f] - g.z[p.i])                                 # Compute direction of travel

        iIdx = getCorrectionIndices(p.i, g.np, g)                       # Initial point correction indices
        fIdx = getCorrectionIndices(p.f, g.np, g)                       # Final point correction indices

        # Choose appropriate branch cuts
        if abs(real(z)) >= abs(imag(z))
            # Up and down moving cuts
            if real(z) >= 0
                θα = -sgn(imag(z)) * π / 2
                θβ =  sgn(imag(z)) * π / 2
            else
                θα = sgn(imag(z)) * π
                θβ = sgn(imag(z)) * π / 2
            end
        else
            # Right and left moving cuts
            if real(z) >= 0
                θα = sgn(imag(z)) * π
                θβ = 0
            else
                θα = 0
                θβ = sgn(imag(z)) * π
            end
        end

        # Trapezoidal weights
        αt = zα.(     g.z[p.p], α, θ = θα)
        βt = zα.(z .- g.z[p.p], β, θ = θβ)

        row[p.p] .+= dir * g.h .* αt .* βt

        # Initial endpoint correction
        h  = n == 1 ? zα(dir, 1 + α, θ = θα) : dir
        αt = n == 1 ? 1 : zα.(g.z[iIdx], α, θ = θα)
        βt = zα.(z .- g.z[iIdx], β, θ = θβ)

        row[iIdx] += h * αt .* βt .* getCorrection(c, dir, n == 1 ? "bp" : "cr")

        # Final endpoint correction
        h  = n == N ? zα(dir, 1 + β, θ = θβ) : dir
        αt = zα.(g.z[fIdx], α, θ = θα)
        βt = n == N ? 1 : zα.(z .- g.z[fIdx], β, θ = θβ)

        row[fIdx] += h * αt .* βt .* getCorrection(c, -dir, n == N ? "ep" : "cr")
    end

    return row
end

"""
    getDiffMat(n, r, α, β; ir = 0.5, er = 5)

Generate differentiation matrix for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` over a grid of radius `r`.

The radius to use the Taylor expansion is given by `ir` while the relative radius of the correction stencils are given by `er`.
"""
function getDiffMat(n, r; α = 0.0, β = 0.0, ir = 0.5, np = 3, nl = 1)
    g = getGrid(n, r, ir = ir, np = np, nl = nl)                        # Generate grid

    (iMap, eMap) = getReducedGridMap(g)                                 # Index maps to reduced grid for indexing through diff matrix

    D = zeros(ComplexF64, length(g.i) + length(g.e), length(g.z))       # Initialize differentiation matrix

    # Populate internal weights using Taylor expansion approximation
    idx = getCorrectionIndices(g.c, g.T, g)                             # Taylor expansion indices
    A = lu(getVand(idx, g))

    for (i, iIdx) ∈ pairs(g.i)
        D[iMap[i], idx] = getInternalWeights(iIdx, A, g, α, β)
    end

    # Populate external weights using generalized Gregory quadrature
    c = getCorrections(g, α = α, β = β)                                 # End/corner corrections

    for (e, eIdx) ∈ pairs(g.e)
        D[eMap[e], :] = getExternalWeights(eIdx, c, g, α, β)
    end

    return D
end
