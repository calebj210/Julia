#=
# Constructors and routines for working with complex grids for quadrature.
#
# Author: Caleb Jacobs
# DLM: September 12, 2023
=#

"Complex grid for use with grid based quadratures."
struct Grid
    "Grid points"
    z::Vector{ComplexF64}

    "Internal indices"
    i::Vector{Int64}

    "External indices"
    e::Vector{Int64}

    "Padding indices"
    p::Vector{Int64}

    "Type of point, 'i', 'e', 'p'"
    t::Vector{Char}

    "Index spacing in x"
    dx::Int             

    "Index spacing in y"
    dy::Int            

    "Origin index"
    c::Int              

    "Grid spacing"
    h::Float64           

    "Radius of square grid"
    r::Float64            
end

"""
    getGrid(n, r; ir = 0.5, p = 0)

Generate a complex grid of radius `r` with `n` nodes from the origin to the adjacent boundaries.

`ir` specifies the radius of the internal nodes and `p` gives the number of padding nodes to add.
"""
function getGrid(n, r; ir = 0.5, p = 0)::Grid
    x⃗ = [range(0, r, length = n + 1)...]                # real parts
    y⃗ = [range(0, r, length = n + 1)...]                # imaginary parts
    h  = abs(x⃗[2] - x⃗[1])                               # Grid spacing

    pad = r .+ h * [1:p...]                             # Padded nodes
    xp = [-pad[end : -1 : 1]; -x⃗[end: -1: 2]; x⃗; pad]   # Padded x vector
    yp = [-pad[end : -1 : 1]; -y⃗[end: -1: 2]; y⃗; pad]   # Padded y vector

    grid = xp' .+ im * yp                               # Padded grid matrix

    dx = stride(grid, 2)                                # Index distance to move in x
    dy = stride(grid, 1)                                # Index distance to move in y
    c  = 1 + (p + n) * (dx + dy)                        # Index of origin

    z⃗ = vec(grid)                                       # Vectorize matrix grid
    i = Array{Int64}([])                                # Initialize internal index array
    e = Array{Int64}([])                                # Initialize external index array
    p = Array{Int64}([])                                # Initialize padding index array
    t = similar(z⃗, Char)                                # Initialize type array

    # Populate index arrays
    for (idx, z) ∈ pairs(z⃗)
        if abs(z) <= ir
            push!(i, idx)                               # Internal node
            t[idx] = 'i'
        elseif abs(real(z)) > r || abs(imag(z)) > r
            push!(p, idx)                               # Padding node
            t[idx] = 'p'
        else
            push!(e, idx)                               # External node
            t[idx] = 'e'
        end
    end

    return Grid(z⃗, i, e, p, t, dx, dy, c, h, r)
end

"""
    getReducedGridMap(g)

Compute index map from internal and external node indices to the non-padded grid.
"""
function getReducedGridMap(g)
    iMap = similar(g.i)             # Interior index map
    eMap = similar(g.e)             # External index map

    iIdx = 1                        # Current internal node index
    eIdx = 1                        # Current external node index
    tIdx = 1                        # Non padded total index

    # Construct the index map
    while iIdx <= length(g.i) && eIdx < length(g.e)
        if g.i[iIdx] < g.e[eIdx]
            iMap[iIdx] = tIdx       # Internal index is next in line
            iIdx += 1               # Move to next internal index
        else
            eMap[eIdx] = tIdx       # External index is next in line
            eIdx += 1               # Move to next external index
        end

        tIdx += 1                   # Move to next index in non-padded grid
    end

    # Populate remaining index maps 
    if iIdx <= length(g.i)
        iMap[iIdx : end] = [tIdx : length(g.i) + length(g.e)...]
    else
        eMap[eIdx : end] = [tIdx : length(g.i) + length(g.e)...]
    end
    
    return(iMap, eMap)
end

"""
    getIdx(z, g)

Find index in `g` of `z`.
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

Get indices to traverse the complex grid from 0 to z with index `zIdx`. The indices are returned as a tuple (horizontal, vertical)
"""
function getPathIndices(zIdx::Int64, g::Grid)
    Nx = round(Int64, real(g.z[zIdx]) / g.h)                        # Number of horizontal nodes
    Ny = round(Int64, imag(g.z[zIdx]) / g.h)                        # Number of vertical nodes
    
    if Nx != 0
#         horiz = g.c .+ 
#                 sign(Nx) * g.dx * [1 : abs(Nx) - 1...]              # March horizontally to z ignoring endpoint
        horiz = g.c + 
                Ny * g.dy .+ 
                sign(Nx) * g.dx * [1 : abs(Nx) - 1...]              # March horizontally to z ignoring endpoint
    else
        horiz = []                                                  # No horizontal nodes
    end

    if Ny != 0
#         vert = g.c + 
#                Nx * g.dx .+
#                sign(Ny) * g.dy * [1 : abs(Ny) - 1...]               # March vertically to z ignoring endpoints
        vert = g.c .+ 
               sign(Ny) * g.dy * [1 : abs(Ny) - 1...]               # March vertically to z ignoring endpoints
    else
        vert = []                                                   # No vertical nodes
    end

    corner = g.c + Ny * g.dy                                        # Index of corner if up then right
#     corner = g.c + Nx * g.dx                                        # Index of corner if right then up

    return (horiz, vert, corner)
end

function getPathIndices(z::ComplexF64, g::Grid)
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
Rotate square correction stencil to point in direction of dir (clockwise)
"""
function rotCorrection(idx::Vector{Int64}, dir::Int64)
    r = length(idx) / 8

    return circshift(idx, -2r * (dir % 4))
end
