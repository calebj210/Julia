#=
# Crank-Nicolson scheme for solving
#   ∂t(u) = ∂x(a(x) ∂x(u)) + f(x,t)
#
# Author: Caleb Jacobs
# DLM: 12-04-2022
=#

using ForwardDiff
using LinearAlgebra
using Plots

default(xlims = (0, 1), ylims = (-1, 1))

"""
    getCRMat(a, x, ht)

Construct Crank-Nicolson matrix given function `a(x)`, uniform grid
data `x`, and time step size `ht`.
"""
function getCRMat(a, x, ht)
    hx = x[2] - x[1]                        # Spatial step size
    
    ax  = a.(x)                             # a(x)  evaluated at x
    adx = ForwardDiff.derivative.(a, x)     # a'(x) evaluated at x

    dl = ht * (-ax[2:end] / (2hx^2) 
       + adx[2:end] / (4hx))                # Lower diagonal

    d  = 1 .+ ht * ax / hx^2                # Diagonal

    du = ht * (-ax[1:end - 1] / (2hx^2) 
       - adx[1:end - 1] / (4hx))            # Upper diagonal

    A = Tridiagonal(dl, d, du)              # Tridiagonal Crank-Nicolson matrix
end

"""
    getRHS(a, f, x, u, t, ht)

Construct right hand side of Crank-Nicolson scheme given functions `a(x)` and 
`f(x,t)`, and data (`x`,`u`) at time `t` with time step size `ht`. 
"""
function getRHS(a, f, x, u, t, ht)
    hx = x[2] - x[1]                        # Spatial step size
    
    ax  = a.(x)                             # a(x)  evaluated at x
    adx = ForwardDiff.derivative.(a, x)     # a'(x) evaluated at x
    
    l = [0; ht * (ax[2:end] / (2hx^2) - adx[2:end] / (4hx)) .* u[1:end - 1]]          # Left node contribution
    m = (1 .- ht * ax / (hx^2)) .* u                                # Center node contribution
    r = [ht * (ax[1:end - 1] / (2hx^2) + adx[1:end - 1] / (4hx)) .* u[2:end]; 0]   # Right node contribution

    return l + m + r + ht * (f.(x, t) + f.(x, t + ht)) / 2
end

"""
    driver(a, f, u0, hx, ht, tf)

Run Crank-Nicolson method for solving model problem.
"""
function driver(a, f, u0, hx, ht, tf)
    x = range(0 + hx, 1 - hx, step = hx)  # Spatial nodes
    X = [0; x; 1]
    t = 0                       # Initialize time
    u = u0.(x)                  # Initialize solution

    A = getCRMat(a, x, ht)
    F(t, u) = getRHS(a, f, x, u, t, ht)
    
    display(plot(X, [0;u;0]))

    while t < tf
        u = A \ F(t, u)
        
        display(plot(X, [0;u;0]))

        t += ht
    end

    display(plot(X, [0;u;0]))

    return u
end
