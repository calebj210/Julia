#=
# FFT Based Solver for Poisson's equation
#   Δu = f
#
# Author: Caleb Jacobs
# DLM: 17-04-2022
=#

using FFTW
using Plots


"""
    getGrid(Nx, Ny)

Generate 2D node grid with `Nx` x-coordinates and `Ny` y-coordinates.
"""
function getGrid(Nx, Ny)
    x = range(0, 1, length = Nx)
    y = range(0, 1, length = Ny)
    
    X = repeat([x...]', Ny, 1)
    Y = repeat(y, 1, Nx)

    return (x, y, X, Y)
end

"""
    poissonSolve(f, Nx, Ny)

Solve Poissons equation over the unit square with homogeneous
boundary Dirichlet boundary condtions and forcing term `f`(x,y).
"""
function poissonSolve(f, Nx, Ny)
    (x, y, X, Y) = getGrid(Nx, Ny)  # Get computation nodes

    fx = map(f, X, Y)               # Function values at each node

    b = FFTW.r2r(fx, FFTW.RODFT00)  # Perform FFT based sine transform

    a = [-b[j, i] / (π^2 * (i^2 + j^2)) 
        for j ∈ 1:Ny, i ∈ 1:Nx]     # Compute FFT of solution

    U = FFTW.r2r(a, FFTW.RODFT00) /
        (4 * (Nx + 1) * (Ny + 1))   # Get solution via inverse sine FFT

    display(contour(x, y, U, st = :contour, 
                    fill = true, 
                    aspect_ratio = :equal))

    return U
end
