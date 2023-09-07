#=
# Generalized Gregory quadrature for computing hypergeometric pFq on a grid
#
# Author: Caleb Jacobs
# DLM: September 7, 2023
=#

using Plots

struct Grid
    z::Vector{ComplexF64}   # Grid points
    dx::Int                 # Index spacing in x
    dy::Int                 # Index spacing in y
    c::Int                  # Index of origin
    h::Float64              # Grid spacing
    r::Float64              # Boundary of grid
end

"""
    plotPath(idx, g)
Display the path defined by idx through the grid g.
"""
function plotPath(idx::Vector{Int64}, g::Grid)
    plt = plot(g.z, st = :scatter, c = :red, ratio = 1, legend = false)
    
    plot!(g.z[idx], st = :scatter, mz = 1 : length(idx), c = :viridis)

    return plt
end
function plotPath(idx::Tuple{Vector, Vector}, g::Grid)
    plt = plot(g.z, st = :scatter, c = :red, ratio = 1, legend = false)
    
    for i ∈ idx
        if !isempty(i)
            plot!(g.z[i], st = :scatter, mz = 1 : length(i), c = :viridis)
        end
    end

    return plt
end

"""
    getGrid(n, r)
Generate a complex grid of radius `r` with n nodes from the origin to the adjacent boundaries.
"""
function getGrid(n, r)
    x⃗ = [range(0, r, length = n + 1)...]                # real parts
    y⃗ = [range(0, r, length = n + 1)...]                # imaginary parts

    grid = [-x⃗[end:-1:2];x⃗]' .+ im * [-y⃗[end:-1:2];y⃗]   # Grid matrix

    dx = stride(grid, 2)                                # Index distance to move in x
    dy = stride(grid, 1)                                # Index distance to move in y
    c  = 1 + n * (dx + dy)                              # Index of origin
    h  = abs(grid[2] - grid[1])                         # Grid spacing

    return Grid(grid[:], dx, dy, c, h, r)
end

"""
    getIdx(z, g)
Find index in g of z.
"""
function getIdx(z::ComplexF64, g::Grid)
    if abs(real(z)) > g.r || abs(imag(z)) > g.r
        throw(DomainError(z, "argument must be inside of the grid"))
    end

    Nx = round(Int64, real(z) / g.h)
    Ny = round(Int64, imag(z) / g.h)

    zIdx = g.c + Nx * g.dx + Ny * g.dy

    if abs(z - g.z[zIdx]) ≈ 0
        return return zIdx
    else
        throw(DomainError(z, "argument must be on the complex grid"))
    end
end

""" 
    getLinearIndices(zIdx, g)
Get indices to traverse the complex grid from 0 to z with index zIdx. The indices are returned as a tuple (horizontal, vertical)
"""
function getLinearIndices(zIdx::Int64, g::Grid)
    Nx = round(Int64, real(g.z[zIdx]) / g.h)
    Ny = round(Int64, imag(g.z[zIdx]) / g.h)
    
    if Nx != 0
        horiz = g.c .+ g.dx * [0 : sign(Nx) : Nx...]            # March horizontally to z
    else
        horiz = []
    end

    if Ny != 0
        vert = g.c + Nx * g.dx .+ g.dy * [0 : sign(Ny) : Ny...] # March vertically to z
    else
        vert = []
    end

    return (horiz, vert)
end

function getLinearIndices(z::ComplexF64, g::Grid)
    return getLinearIndices(getIdx(z, g), g)
end

""" 
    getCorrectionIndices(z, r, z⃗)
Get indices of surrounding correction nodes in z⃗ an adjacent distance of r away from z.
"""
function getCorrectionIndices(zIdx::Int64, r::Int64, g::Grid)
    if r == 0
        return [zIdx]
    end

    idx = zeros(Int64, 8r)                                              # Initialize indices

    idx[     1 :  r + 1] = zIdx + r * g.dx .+ [     0 : r...] * g.dy    # Right top wall
    idx[ r + 2 : 3r + 1] = zIdx + r * g.dy .- [-r + 1 : r...] * g.dx    # Top wall
    idx[3r + 2 : 5r + 1] = zIdx - r * g.dx .- [-r + 1 : r...] * g.dy    # Left wall
    idx[5r + 2 : 7r + 1] = zIdx - r * g.dy .+ [-r + 1 : r...] * g.dx    # Bottom wall

    if r > 1
        idx[7r + 2 : 8r] = zIdx + r * g.dx .+ [-r + 1 : - 1...] * g.dy  # Right bottom wall
    end

    return idx
end

function getCorrectionIndices(z::ComplexF64, r::Int64, g::Grid)
    return getCorrectionIndices(getIdx(z, g), r, g)
end

"""
    rotCorrection(idx, dir, g::Grid)
Rotate square correction stencil to point in direction of dir
"""
function rotCorrection(idx::Vector{Int64}, dir::Int64, g::Grid)
    r = length(idx) / 8

    return circshift(idx, -2r * (dir % 4))
end
