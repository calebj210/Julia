#=
#   ODE approach to computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: June 6, 2024
=#

using ForwardDiff: derivative
include("TimeStep.jl")

compDiff(f, z, n = 1) = 
    derivative(x -> real(f(z + x)), 0) + derivative(x -> imag(f(z + x)), 0) * im

function taylorA(a, b, c, z, tol = 1e-15)
    if c == round(c, RoundToZero) && real(c) <= 0
        a1 = gamma(a - c + 1) / gamma(a) * gamma(b - c + 1) / gamma(b) / gamma(2 - c) * z^(1 - c)
        b1 = a1

        for j ∈ 1 : 500
            a1 = (a - c + j) * (b - c + j) / j * z / (j - c + 1) * a1
            b1 += a1
        end
    else
        a1 = b1 = 1

        for j = 1 : 500
            a1 = (a + j - 1) * (b + j - 1) / (c + j - 1) * z / j * a1
            b1 += a1
        end
    end


    return b1
end

function F21(a, b, c, z::Number; z0 = .25 + 0im, N = 10)
    fp = compDiff(x -> taylorA(a, b, c, x), z0)
    f0 = [taylorA(a, b, c, z), fp]
    display(f0)

    F(zn, fn) = [fn[2], (a * b * fn[1] - (c - (a + b + 1) * zn) * fn[2]) / zn / (1 - zn)]

    f = odeSolve(z0, f0, F, z, N = N)[1]

    return f
end
