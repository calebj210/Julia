#=
#   ODE approach to computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: June 6, 2024
=#

include("TimeStep.jl")

function F21(a, b, c, z::Number; z0 = .25 + 0im, f0 = [1.0 + 0im, .1 + 0im], N = 10)
    F(zn, fn) = [fn[2], (a * b * fn[1] - (c - (a + b + 1) * zn) * fn[2]) / zn / (1 - zn)]

    f = odeSolve(z0, f0, F, z, N = N)[1]

    return f
end
