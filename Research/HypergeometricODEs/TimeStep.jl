#=
# Time stepping algorithms for complex steps in ODEs
#
# Author: Caleb Jacobs
# DLM: June 20, 2024
=#

using Polynomials

"""
    /(p::AbstractPolynomial, q::AbstractPolymomial)
Compute the polynomial quotient p / q tossing the remainder away.
"""
function Base.:/(a::AbstractPolynomial, b::AbstractPolynomial)
    n = length(a)                   # Number of coefficients to compute
    x = variable(a)                 # Base variable

    q, r = divrem(a, b)             # Quotient and remainder
    
    d = Polynomial(reverse(b.coeffs))
    r = Polynomial(reverse(r.coeffs)) * x^(length(d) - length(r))

    # Expand remainder and add it to quotient
    for i ∈ 0 : n
        quotient, r = divrem(r, d)
        r *= x
        q += quotient * x^(i)
    end

    return q
end

function mySetindex!(p::AbstractPolynomial, value, idx::Int)
    n = length(coeffs(p))
    if n ≤ idx
        resize!(p.coeffs, idx + 1)
        p.coeffs[n + 1:idx] .= 0
    end
    p.coeffs[idx + 1] = value
    return p
end

function rk4Step(f, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h / 2, y + h * k1 / 2)
    k3 = f(t + h / 2, y + h * k2 / 2)
    k4 = f(t + h, y + h * k3)

    return h * (k1 + 2k2 + 2k3 + k4) / 6
end

"""
    taylorStep(z0, y0, f, h; order = 10)

Compute Taylor expansion ODE step for the IVP y'(z) = `f`(z, y), y(`z0`) = `y0` given a step size of `h`.
"""
function taylorStep(z0, y0, f, h; order = 10)
    z = Polynomial([z0, 1])                         # Time step polynomial
    y = Polynomial.(y0)                             # Solution polynomial

    for i ∈ 1 : order
        RHS = integrate.(f(z, y))                   # Compute next order of expansion

        for j ∈ eachindex(y)
            mySetindex!(y[j], RHS[j][i], i)         # Store next expansion coefficients
        end
    end

    return [p(h) for p ∈ y]                         # Evaluate expansion
end

"""
    odeSolve(z0, f0, f, z, H; order = 20)

Solve IVP ODE y'(z) = `f`(z, y) with initial conditions y(`z0`) = `f0` using a max step size of `H` and taylor expansions up to degree of `order`.
"""
function odeSolve(z0, y0, f, z, H; order = 20, store = true)
    dir = sign(z - z0)                                  # Step size

    tn = z0                                             # Time variable
    yn = copy(y0)                                       # ODE solution

    if store
        Z  = [tn]                                       # Time vector
        Y  = [yn]                                       # Solution vector
    end

    for i ∈ 1 : ceil(Int64, abs(z - z0) / H)
        h = dir * min(H, abs(z - tn))                   # Min step size allowed
        yn = taylorStep(tn, yn, f, h, order = order)    # March with Taylor series
        tn += h                                         # Step time     

        if store
            push!(Z, tn)                                # Next z value
            push!(Y, yn)                                # Next solution value
        end
    end

    if store
        return (Z, Y) 
    else
        return yn
    end
end
