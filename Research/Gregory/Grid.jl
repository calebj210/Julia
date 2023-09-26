#=
# Constructors and routines for working with complex grids for quadrature.
#
# Author: Caleb Jacobs
# DLM: September 25, 2023
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

    "Number of padded nodes"
    np::Int64
    
    "Index radius of Taylor expansion"
    T::Int64
end

struct Path
    i::Int64
    p::Vector{Int64}
    f::Int64
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

    np = p                                              # Number of padding nodes

    pad = r .+ h * [1:p...]                             # Padded nodes
    xp = [-pad[end : -1 : 1]; -x⃗[end: -1: 2]; x⃗; pad]   # Padded x vector
    yp = [-pad[end : -1 : 1]; -y⃗[end: -1: 2]; y⃗; pad]   # Padded y vector

    grid = xp' .+ im * yp                               # Padded grid matrix

    dx = stride(grid, 2)                                # Index distance to move in x
    dy = stride(grid, 1)                                # Index distance to move in y
    c  = 1 + (p + n) * (dx + dy)                        # Index of origin
    T = round(Int64, ir / h)

    z⃗ = vec(grid)                                       # Vectorize matrix grid
    i = Array{Int64}([])                                # Initialize internal index array
    e = Array{Int64}([])                                # Initialize external index array
    p = Array{Int64}([])                                # Initialize padding index array

    # Populate index arrays
    for (idx, z) ∈ pairs(z⃗)
        if abs(z) <= ir
            push!(i, idx)                               # Internal node
        elseif abs(real(z)) > r || abs(imag(z)) > r
            push!(p, idx)                               # Padding node
        else
            push!(e, idx)                               # External node
        end
    end

    return Grid(z⃗, i, e, p, dx, dy, c, h, r, np, T)
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
    getPath(zIdx, g)

Get indices to traverse the complex grid from 0 to z with index `zIdx`. The indices are returned as a tuple (horizontal, vertical)
"""
function getPath(zIdx::Int64, g::Grid, r)
    Nx = round(Int64, real(g.z[zIdx]) / g.h)                        # Number of horizontal nodes
    Ny = round(Int64, imag(g.z[zIdx]) / g.h)                        # Number of vertical nodes

    sgn(x) = x >= 0 ? 1 : -1                                        # Nonzero returning sign function

    # Begin u-pathing
    if abs(Nx) >= abs(Ny)
        v1 = g.c .- g.dy * sgn(Ny) * [1 : 2r - 1...]

        c1 = v1[end] - g.dy * sgn(Ny)

        h  = c1 .+ g.dx * sgn(Nx) * [1 : abs(Nx) - 1...]

        c2 = h[end] + g.dx * sgn(Nx)

        v2 = c2 .+ g.dy * sgn(Ny) * [1 : 2r + abs(Ny) - 1...]

        p1 = Path(g.c, v1, c1)
        p2 = Path(c1,  h,  c2)
        p3 = Path(c2,  v2, zIdx)

        return [p1, p2, p3]
    else
        h1 = g.c .- g.dx * sgn(Nx) * [1 : 2r - 1...]

        c1 = h1[end] - g.dx * sgn(Nx)

        v  = c1 .+ g.dy * sgn(Ny) * [1 : abs(Ny) - 1...]

        c2 = v[end] + g.dy * sgn(Ny)

        h2 = c2 .+ g.dx * sgn(Nx) * [1 : 2r + abs(Nx) - 1...]

        p1 = Path(g.c, h1, c1)
        p2 = Path(c1,  v,  c2)
        p3 = Path(c2,  h2, zIdx)

        return [p1, p2, p3]
    end

    return path
end

function getPathIndices(z::ComplexF64, g::Grid)
    return getLinearIndices(getIdx(z, g), g)
end

""" 
    getCorrectionIndices(z, r, g)
Get indices of surrounding correction nodes in grid g an adjacent distance of r away from z.
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
function rotCorrection(idx::Vector, dir::Int64)
    r = length(idx) / 8

    return circshift(idx, -2r * (dir % 4))
end
