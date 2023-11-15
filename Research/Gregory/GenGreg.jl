#=
# Generalized, grid-based Gregory quadrature for computing hypergeometric pFq
#
# Author: Caleb Jacobs
# DLM: November 9, 2023
=#

using SpecialFunctions
using SparseArrays
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

"Modified sign function to return 1 when z = 0"
sgn(z) = iszero(z) ? one(z) : sign(z)

"Modified sign function to return -1 when z = 0"
nsgn(z) = iszero(z) ? -one(z) : sign(z)

"Compute roots given a power α and a branch cut rotation of θ."
function zα(z, α::Real; θ = 0)
    if θ >= 0
        return angle(z) >= θ - π ? z^α : z^α * cispi( 2(α % 1))
    else
        return angle(z) <  θ + π ? z^α : z^α * cispi(-2(α % 1))
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

"Get extra branch rotation correction"
function θγ(z, ze, γ::Real)
    if imag(z) == 0
        return imag(ze) >= 0 ? cispi(2(γ % 1)) : 1
    elseif sgn(imag(z)) == sgn(imag(ze))
        return 1
    else
        return cispi(-2sgn(imag(z)) * (γ % 1))
    end
end

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
    getBranchAngles(z)

Compute appropriate branch cut rotation angles based on the path to `z`.
"""
function getBranchAngle(z, g::Grid; branch = false)
    if !branch
        if abs(real(z)) < abs(imag(z)) || real(z) >= 0 && abs(imag(z)) >= 3g.np * g.h
            # Right and left moving cuts
            if real(z) >= 0
                θα = sgn(imag(z)) * π
                θβ = 0
            else
                θα = 0
                θβ = sgn(imag(z)) * π
            end
        else
            # Up and down moving cuts
            if real(z) >= 0
                if iszero(imag(z))
                    θα = 0 
                    θβ = π / 2
                else
                    θα = 0
                    θβ = -sgn(imag(z)) * π / 2
                end
            else
                θα = sgn(imag(z)) * π
                θβ = sgn(imag(z)) * π / 2
            end
        end
    else
        if iszero(imag(z))
            θα = 0 
            θβ = -π / 2
        else
            θα = 0
            θβ = nsgn(imag(z)) * π / 2
        end
    end

    return (θα, θβ)
end

# """
#     getExternalWeights(row, zIdx, path, g::Grid, α, β)
# 
# Get external differentiation entries corresponding to `g`.z[`zIdx`].
# """
# function getExternalWeights(zIdx, c::Corrections, g::Grid, α, β)
#     path = getPath(zIdx, g, 2g.np)                                      # Get path from origin to node
#     N = length(path)                                                    # Number of paths to travel
#     z = g.z[zIdx]                                                       # Current z value
#     row = zeros(ComplexF64, length(g.z))                                # Initialize row to populate
# 
#     for (n, p) ∈ pairs(path)
#         dir = sign(g.z[p.f] - g.z[p.i])                                 # Compute direction of travel
# 
#         iIdx = getCorrectionIndices(p.i, g.np, g)                       # Initial point correction indices
#         fIdx = getCorrectionIndices(p.f, g.np, g)                       # Final point correction indices
# 
#         # Choose appropriate branch cuts
#         θα, θβ = getBranchAngle(z, g)
# 
#         # Trapezoidal weights
#         αt = zα.(     g.z[p.p], α, θ = θα)
#         βt = zα.(z .- g.z[p.p], β, θ = θβ)
# 
#         row[p.p] .+= dir * g.h .* αt .* βt
# 
#         # Initial endpoint correction
#         h  = n == 1 ? zα(dir, 1 + α, θ = θα) : dir
#         αt = n == 1 ? 1 : zα.(g.z[iIdx], α, θ = θα)
#         βt = zα.(z .- g.z[iIdx], β, θ = θβ)
# 
#         row[iIdx] += h * αt .* βt .* getCorrection(c, dir, n == 1 ? "bp" : "cr")
# 
#         # Final endpoint correction
#         h  = n == N ? zα(dir, 1 + β, θ = θβ) : dir
#         αt = zα.(g.z[fIdx], α, θ = θα)
#         βt = n == N ? 1 : zα.(z .- g.z[fIdx], β, θ = θβ)
# 
#         row[fIdx] += h * αt .* βt .* getCorrection(c, -dir, n == N ? "ep" : "cr")
#     end
#     
#     tmp = .!iszero.(row)
#     rowIdx = findall(tmp)
#     row = row[tmp]
# 
#     return (rowIdx, row)
# end

function sortPath(idx, g, z; branch = false)
    pIdx = Vector{Int64}()
    bIdx = Vector{Int64}()

    # No need for branching
    if real(z) < 1
        return ([1:length(idx)...], bIdx)
    end

    # Decide whether the index is one sheet or another
    for (i, v) ∈ pairs(idx)
        if !branch
            if nsgn(imag(g.z[v])) == nsgn(imag(z))
                push!(pIdx, i)
            else
                push!(bIdx, i)
            end
        else
            if imag(z) == 0
                if nsgn(imag(g.z[v])) == 1
                    push!(pIdx, i)
                else
                    push!(bIdx, i)
                end
            elseif imag(z) > 0
                if nsgn(imag(g.z[v])) == -1
                    push!(pIdx, i)
                else
                    push!(bIdx, i)
                end
            else
                if nsgn(imag(g.z[v])) == 1
                    push!(pIdx, i)
                else
                    push!(bIdx, i)
                end
            end
        end

#         if (nsgn(imag(g.z[v])) == nsgn(imag(z))) ⊻ branch
#             push!(pIdx, i)
#         else
#             push!(bIdx, i)
#         end
    end

    return (pIdx, bIdx)
end

function getExternalWeights(zIdx, c::Corrections, g::Grid, α, β; branch = false)
    path = getPath(zIdx, g, 2g.np, branch = branch)                     # Get path from origin to node
    N = length(path)                                                    # Number of paths to travel
    z = g.z[zIdx]                                                       # Current z value

    rowf = zeros(ComplexF64, length(g.z))                               # Initialize main row to populate
    rowh = zeros(ComplexF64, length(g.z))                               # Initialize branch row to populate

    for (n, p) ∈ pairs(path)
        dir = sign(g.z[p.f] - g.z[p.i])                                 # Compute direction of travel

        iIdx = getCorrectionIndices(p.i, g.np, g)                       # Initial point correction indices
        fIdx = getCorrectionIndices(p.f, g.np, g)                       # Final point correction indices

        # Fix indices if needed
        if n == N && real(z) > 1
            ppIdx, bpIdx = sortPath(p.p,  g, z, branch = branch)
            piIdx, biIdx = sortPath(iIdx, g, z, branch = branch)
            pfIdx, bfIdx = sortPath(fIdx, g, z, branch = branch)
        else
            ppIdx, bpIdx = [1 : length(p.p)...],  []
            piIdx, biIdx = [1 : length(iIdx)...], []
            pfIdx, bfIdx = [1 : length(fIdx)...], []
        end

        # Choose appropriate branch cuts
        θα, θβ = getBranchAngle(z, g, branch = branch)

        # Trapezoidal weights
        αt = zα.(     g.z[p.p], α, θ = θα)
        βt = zα.(z .- g.z[p.p], β, θ = θβ)

        rowf[p.p[ppIdx]] .+= dir * g.h .* αt[ppIdx] .* βt[ppIdx]
        rowh[p.p[bpIdx]] .+= dir * g.h .* αt[bpIdx] .* βt[bpIdx]

        # Initial endpoint correction
        h  = n == 1 ? zα(dir, 1 + α, θ = θα) : dir
        αt = n == 1 ? ones(length(iIdx)) : zα.(g.z[iIdx], α, θ = θα)
        βt = zα.(z .- g.z[iIdx], β, θ = θβ)
        
        tmp = getCorrection(c, dir, n == 1 ? "bp" : "cr")

        rowf[iIdx[piIdx]] += h * αt[piIdx] .* βt[piIdx] .* tmp[piIdx]
        rowh[iIdx[biIdx]] += h * αt[biIdx] .* βt[biIdx] .* tmp[biIdx]
        
        # Final endpoint correction
        h  = n == N ? zα(dir, 1 + β, θ = θβ) : dir
        αt = zα.(g.z[fIdx], α, θ = θα)
        βt = n == N ? ones(length(fIdx)) : zα.(z .- g.z[fIdx], β, θ = θβ)

        tmp = getCorrection(c, -dir, n == N ? "ep" : "cr")

        rowf[fIdx[pfIdx]] += h * αt[pfIdx] .* βt[pfIdx] .* tmp[pfIdx]
        rowh[fIdx[bfIdx]] += h * αt[bfIdx] .* βt[bfIdx] .* tmp[bfIdx]
    end

    tmp = .!iszero.(rowf)
    rowfIdx = findall(tmp)
    rowf = rowf[tmp]

    tmp = .!iszero.(rowh)
    rowhIdx = findall(tmp)
    rowh = rowh[tmp]

    return (rowfIdx, rowf, rowhIdx, rowh)
end

"""
    getDiffMat(n, r, α, β; ir = 0.5, er = 5)

Generate differentiation matrix for ``∫₀ᶻ(u)ᵅ(z-u)ᵝf(u)du`` over a grid of radius `r`.

The radius to use the Taylor expansion is given by `ir` while the relative radius of the correction stencils are given by `er`.
"""
function getDiffMat(n, r; α = 0.0, β = 0.0, ir = 0.5, np = 3, nl = 1, branch = false)
    g = getGrid(n, r, ir = ir, np = np, nl = nl)                            # Generate grid

#     (iMap, eMap) = getReducedGridMap(g)                                     # Index maps to reduced grid for indexing through diff matrix
    M = length(g.i) + length(g.e)
    N = length(g.z)

    if !branch
        row0 = Vector{Int64}()
        row1 = Vector{Int64}()

        col0 = Vector{Int64}()
        col1 = Vector{Int64}()
        
        wgt0 = Vector{ComplexF64}()
        wgt1 = Vector{ComplexF64}()
    else
        # Initialize row vectors
        row0 = Vector{Int64}()
        row1 = Vector{Int64}()
        row2 = Vector{Int64}()
        row3 = Vector{Int64}()

        # Initialize column vectors
        col0 = Vector{Int64}()
        col1 = Vector{Int64}()
        col2 = Vector{Int64}()
        col3 = Vector{Int64}()

        # Initialize weight vectors
        wgt0 = Vector{ComplexF64}()
        wgt1 = Vector{ComplexF64}()
        wgt2 = Vector{ComplexF64}()
        wgt3 = Vector{ComplexF64}()
    end

    # Populate internal weights using Taylor expansion approximation
    A = lu(getVand(g.ib, g))                                                # Compute vandermonde of internal boundary nodes

    # Populate external weights using generalized Gregory quadrature
    c = getCorrections(g, α = α, β = β)                                     # End/corner corrections

    i = 1
    for (j, t) ∈ pairs(g.ie)
        if t == 'p'
            continue
        end

        if t == 'i'
            push!(row0, g.ib...)
            push!(col0, repeat([i], length(g.ib))...)
            push!(wgt0, getInternalWeights(j, A, g, α, β)...)
            
            i += 1
        elseif t == 'e'
            idx0, wgts0, idx1, wgts1 = getExternalWeights(j, c, g, α, β, branch = false)

            push!(row0, idx0...)
            push!(col0, repeat([i], length(idx0))...)
            push!(wgt0, wgts0...)

            push!(row1, idx1...)
            push!(col1, repeat([i], length(idx1))...)
            push!(wgt1, wgts1...)
            if branch && real(g.z[i]) >= 1
                idx2, wgts2, idx3, wgts3 = getExternalWeights(j, c, g, α, β, branch = true)

                push!(row2, idx2...)
                push!(col2, repeat([i], length(idx2))...)
                push!(wgt2, wgts2...)

                push!(row3, idx3...)
                push!(col3, repeat([i], length(idx3))...)
                push!(wgt3, wgts3...)
            end

            i += 1
        end
    end

    # Construct sparse arrays
    D0 = sparse(col0, row0, wgt0, M, N)
    D1 = sparse(col1, row1, wgt1, M, N)
    if branch
        D2 = sparse(col2, row2, wgt2, M, N)
        D3 = sparse(col3, row3, wgt3, M, N)
    end

    if !branch
        return (D0, D1)
    else
        return (D0, D1, D2, D3)
    end
end
