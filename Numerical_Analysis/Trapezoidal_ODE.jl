#=
# 2×2 First Order ODE Solver Using Trapezoidal rule
#
# Author: Caleb Jacobs
# Date last modified: 14-Mar-2022
=#

using Plots
using SpecialFunctions
using ForwardDiff
using LinearAlgebra

# First order Bessel Equation system
function f(t, y)
    # Use approximation if we are near the singularity t = 0
    if t < 1e-10
        return [y[2], -3*t/8] # Higher order terms + 5*(t^3)/96 - 7*(t^5)/3072]
    else
        return [y[2], ((1 - t^2)*y[1] - t*y[2]) / (t^2)]
    end
end

# Newton method system solver
function newton(f; maxIts = 100, ε = 1e-8, y0 = [0, 0])
    y = y0                  # Initial guess
    
    for i ∈ 1 : maxIts
        J = ForwardDiff.jacobian(f, y)  # Get jacobian
        ynew = y - (J \ f(y))           # Find next iterate

        # Check for convergence
        if norm(y - ynew) <= ε
            y = ynew

            return y
        end

        y = ynew                        # Pass to next iteration
    end

    return y
end

# Trapezoidal rule
function trapz(f, a, y0, h, n)
    yi = y0             # Initial conditions
    ti = a              # Initial time
    tf = a + h          # First time step

    # Run trapezoidal until desired time
    for i = 1 : n
        # Current trapezoidal equation
        g(y) = y - (yi + h * (f(ti, yi) + f(tf, y)) / 2)

        yi = newton(g)  # Solve trapezoidal equation

        ti = tf         # Store new time
        tf = ti + h     # Compute next time
    end

    return yi
end

# Trapezoidal rule with Richardson Extrapolation
function richTrap(f, a, b, y0; n = 1, rn = 1)
    h = (b - a) / n                 # Compute time step
    
    r = zeros(Float64, rn, rn)      # Initialize richardson matrix
    sol = trapz(f, a, y0, h, n)     # Get initial solution
    r[1, 1] = sol[1]                # Store initial solution

    # Begin Richardson exptrapolation
    for i = 1 : rn - 1
        h /= 2                      # Half time step size
        n *= 2                      # Double number of step to take
        
        sol = trapz(f, a, y0, h, n) # Get solution with current step size
        r[i + 1, 1] = sol[1]        # Store solution

        # Compute richardson exptrapolation with current data
        for j = 1 : i
            r[i + 1, j + 1] = ((4^j) * r[i + 1, j] - r[i, j]) / (4^j - 1)
        end
    end

    return r[rn, rn]
end

besselJ = richTrap(f, 0, 3*π, [0,1/2], n = 40, rn = 10)
display(besselJ)
display(besselJ - besselj(1, 3*π))
