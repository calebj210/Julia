#=
# Central Difference Based Solver for
#   ∂tt(u) = ∂x(a(x) ∂x(u)) + f(x,t)
#
# Author: Caleb Jacobs
# DLM: 12-04-2022
=#

using Plots
using ForwardDiff

default(xlims = (0,1), ylims = (-1, 1))

"""
    solveInitial(x, hx, ht, u0, u1)

Solve for first time step incorporating boundary data `u0` and `u1`.
"""
function solveInitial(a, f, x, u0, u1, ht)
    n = size(x, 1)                                      # Number of nodes
    t = 0                                               # Initial time
    hx = x[2] - x[1]                                    # Spatial stepsize

    r = (ht^2 / hx^2) * a.(x)                           # a(x)  evaluated at x
    k = (ht^2 / (2hx)) * ForwardDiff.derivative.(a, x)  # a'(x) evaluated at x
    u1x = u1.(x)                                        # u1(x) evaluated at x

    u = zeros(n, 2)                                     # Initialize solution
    u[:, 1] = u0.(x)                                    # Initial condition

    inr = 2:(n - 1)                                     # Inner range
    otr = [1, n]                                        # Outer range

    display(plot(x, u[:, 1]))
    sleep(1)
    
    # Compute inner node step
    u[inr, 2] .= ((r[inr] + k[inr]) .* u[inr .+ 1, 1] +
                 (2 .- 2r[inr])    .* u[inr, 1]       +
                 (r[inr] - k[inr]) .* u[inr .- 1, 1]  +
                  2ht * u1x[inr] + ht^2 * f.(x[inr], t)) / 2

    # Compute boundary node step using period BCs
    u[otr, 2] .= ((r[otr] + k[otr])  * u[2, 1]    .+
                 (2 .- 2r[otr])    .* u[otr, 1]   .+
                 (r[otr] - k[otr])  * u[n - 1, 1]  +
                  2ht * u1x[otr] + ht^2 * f.(x[otr], t)) / 2

    return u
end

"""
    solveFD(a, f, hx, ht, u0, u1, tf)

Solve wave-like problem given standard constraints.
"""
function solveFD(a, f, hx, ht, u0, u1, tf)
    x = range(0, 1, step = hx)                          # Spatial nodes
    n = size(x, 1)                                      # Number of nodes
    t = 0                                               # Initialize time

    u = solveInitial(a, f, x, u0, u1, ht)               # Initial solution
    uNew = zeros(n)                                     # Initialize solution vector

    r = (ht^2 / hx^2) * a.(x)                           # a(x)  evaluated at x
    k = (ht^2 / (2hx)) * ForwardDiff.derivative.(a, x)  # a'(x) evaluated at x

    inr = 2:(n - 1)                                     # Inner range
    otr = [1, n]                                        # Outer range
    
    display(plot(x, u[:,1]))

    while t < tf
        # Compute inner node step
        uNew[inr] .= (r[inr] + k[inr]) .* u[inr .+ 1, 2] +
                     (2 .- 2r[inr])    .* u[inr, 2]      +
                     (r[inr] - k[inr]) .* u[inr .- 1, 2] -
                     u[inr, 1]  +  ht^2 * f.(x[inr], t)

        # Compute boundary node step using period BCs
        uNew[otr] .= (r[otr] + k[otr])  * u[2, 2]     .+
                     (2 .- 2r[otr])    .* u[otr, 2]   .+
                     (r[otr] - k[otr])  * u[n - 1, 2]  -
                     u[otr, 1]  +  ht^2 * f.(x[otr], t)

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
