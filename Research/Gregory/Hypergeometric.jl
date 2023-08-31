#=
# Routines for computing hypergeometric functions
#
# Author: Caleb Jacobs
# DLM: August 31, 2023
=#

using SpecialFunctions
include("Gregory.jl")

"""
    _2F_1(a, b, c, z)
Compute the ₂F₁(a, b; c; z) hypergeometric function using DLMF 15.6.1
"""
function _2F_1(a, b, c, z)
    α = 1 - b
    β = 1 + b - c

    if z == 0
        return dftInt(t -> 1, 0, 1, α = α, β = β)
    elseif z == 1
        return Inf + im * Inf
    end

    val = dftInt(t -> (1 - z * t)^(-a), 0, 1, α = α, β = β)

    return val / gamma(b) / gamma(c - b)
end
