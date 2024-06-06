#=
# Time stepping algorithms for complex steps in ODEs
#
# Author: Caleb Jacobs
# DLM: June 5, 2024
=#

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

function myExp(z; N = 10)
    h = z / N


    f(x, y) = -y[1] * x
    F(x, y) = [y[2], f(x, y)]

    yn = [0.0 + 0im, 1.0 + 0im]
    tn = 0.0 + 0.0im
    for i ∈ 1 : N
        yn += rk4Step(F, tn, yn, h) 
        tn += h
    end

    return yn[1]
end
