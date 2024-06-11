#=
# Time stepping algorithms for complex steps in ODEs
#
# Author: Caleb Jacobs
# DLM: June 10, 2024
=#

using Polynomials

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

function odeSolve(z0, f0, F, z; N = 10)
    h = (z - z0) / N                    # Step size

    tn = z0
    yn = copy(f0)
    for i ∈ 1 : N
        f0 += rk4Step(F, tn, yn, h)
        tn += h
    end

    return f0
end

function myExp(z; N = 10, order = 10)
    h = z / N

    F(x, y) = [y[2], -y[1]]

    yn = [0.0 + 0im, 1.0 + 0im]
    tn = 0.0 + 0.0im
    for i ∈ 1 : N
        yn = taylor(tn, yn, F, h, order = order) 
        tn += h
    end

    return yn[1]
end

function taylor(z0, y0, f, h; order = 10)
    z = Polynomial([z0, 1], :h)
    y = Polynomial.(y0, :h)

    for i ∈ 1 : order
        RHS = integrate.(f(z, y))

        for j ∈ eachindex(y)
            mySetindex!(y[j], RHS[j][i], i)
        end
    end


    return [p(h) for p ∈ y]
end
