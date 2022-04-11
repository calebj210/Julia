#=
# Central Difference Based Solver for
#   ∂tt(u) = ∂x(a(x) ∂x(u)) + f(x,t)
#
# Author: Caleb Jacobs
# DLM: 10-04-2022
=#

using Plots
using ForwardDiff

"""
    solveInitial(x, hx, ht, u0, u1)

Solve for first time step incorporating boundary data `u0` and `u1`.
"""
function solveInitial(a, f, x, t, ht, u0, u1)
    n = size(x, 1)                      # Number of nodes

    hx = x[2] - x[1]                    # Spatial stepsize

    ax  = a.(x)                         # a(x)  evaluated at x
    adx = ForwardDiff.derivative.(a, x) # a'(x) evaluated at x
    u1x = u1.(x)                        # u1(x) evaluated at x

    u = zeros(n, 2)                     # Initialize solution
    u[:, 1] = u0.(x)                    # Initial condition

    inr = 2:(n - 1)                     # Inner range

    display(plot(x, u[:, 1]))
    sleep(1)
    
    # Compute inner node step
    u[inr, 2] =   (2ht * u1x[inr] 
                + (ax[inr] * ht^2 / (hx^2) + adx[inr] * ht^2 / (2hx)) .* u[inr .+ 1, 1] 
                + (2 .- ax[inr] * 2ht^2 / (hx^2)) .* u[inr, 1] 
                + (ax[inr] * ht^2 / (hx^2) - adx[inr] * ht^2 / (2hx)) .* u[inr .- 1, 1]
                + ht^2 * f.(x[inr], t)) / 2

    # Compute boundary node step using period BCs
    u[[1,n], 2] .=   (2ht * u1x[1]
                   + (ax[1] * ht^2 / (hx^2) + adx[1] * ht^2 / (2hx)) * u[2, 1] 
                   + (2 - ax[1] * 2ht^2 / (hx^2)) * u[1, 1] 
                   + (ax[1] * ht^2 / (hx^2) - adx[1] * ht^2 / (2hx)) * u[n - 1, 1]
                   + ht^2 * f.(x[1], t)) / 2

    return u
end

"""
    solveFD(a, f, hx, ht, u0, u1, tf)

Solve wave-like problem given standard constraints.
"""
function solveFD(a, f, hx, ht, u0, u1, tf)
    x = range(0, 1, step = hx)          # Spatial nodes
    n = size(x, 1)                      # Number of nodes
    t = 0                               # Initialize time

    u = solveInitial(a, f, x, t, ht, u0, u1) # Initial solution
    uNew = zeros(n)                     # Initialize solution vector

    ax  = a.(x)                         # a(x)  evaluated at x
    adx = ForwardDiff.derivative.(a, x) # a'(x) evaluated at x

    inr = 2:(n - 1)                     # Inner range
    
    display(plot(x, u[:,1]))

    while t < tf
        # Compute inner node step
        uNew[inr] .=  ((ax[inr] * ht^2 / (hx^2) + adx[inr] * ht^2 / (2hx)) .* u[inr .+ 1, 2] 
                    + (2 .- ax[inr] * 2ht^2 / (hx^2)) .* u[inr, 2] 
                    + (ax[inr] * ht^2 / (hx^2) - adx[inr] * ht^2 / (2hx)) .* u[inr .- 1, 2]
                    + ht^2 * f.(x[inr], t)
                    - u[inr, 2])
    
        # Compute boundary node step using period BCs
        uNew[[1,n]] .=   ((ax[1] * ht^2 / (hx^2) + adx[1] * ht^2 / (2hx)) * u[2, 2] 
                       + (2 - ax[1] * 2ht^2 / (hx^2)) * u[1, 2] 
                       + (ax[1] * ht^2 / (hx^2) - adx[1] * ht^2 / (2hx)) * u[n - 1, 2]
                       + ht^2 * f.(x[1], t)
                       - u[1, 1])

        u[:, 1] = u[:, 2]   # Move current nodes back
        u[:, 2] = uNew      # Move new nodes into current
        t += ht             # Update time

        display(plot(x, uNew))
    end

    return uNew
end

function driver(a, f, hx, ht, u0, u1, tf)
    sol = solveFD(a, f, hx, ht, u0, u1, tf)
end
