#=
# Simple gregory integration
#
=#

using SpecialFunctions
using LinearAlgebra

"""
    vand(z⃗)
Construct a rotated vandermonde matrix using the elements in `z⃗`.
"""
function vand(z⃗ :: Vector)
    V = [z^i for i ∈ 0 : length(z⃗) - 1, z ∈ z⃗]

    return V
end

"""
    ζ⃗(α, h, z⃗)
Construct RHS of general trapezoidal rule endpoint correction.
"""
function ζ⃗(α, h, z⃗)
    RHS = [zeta(α + 1 - i) * h^i for i ∈ 0 : length(z⃗) - 1]

    return RHS
end

"""
    getCorrections(α, h, z⃗)
Compute trapezoidal rule end corrections given `z⃗` is a list of the end correction nodes relative to the base point, `α` is the order of the singularity
"""
function getCorrections(α, h, z⃗)
    w = vand(z⃗) \ ζ⃗(α, h, z⃗)

    return w
end

"""
    nonSingWeights(h, z⃗)
Compute singularity free end corrections given spacing `h` and correction nodes `z⃗`.
"""
function nonSingWeights(h, z⃗)
    w = -getCorrections(-1, h, z⃗)   # Generate most of the weights
    w[1] -= 1                       # Fix first weight

    return w
end

"""
    nonSingInt(f, a, b, n, p)
Integrate a function `f` from `a` to `b` with no singularities using `n` nodes and order `p`
"""
function nonSingInt(f, a, b, n, p)
    z⃗ = [range(a, b, length = n)...]                # Trapezoidal nodes
    h = z⃗[2] - z⃗[1]                                 # Step size

    w = ones(n)                                     # Uncorrected TR weights
    if p >= 2                                       # Compute corrections if needed
        c = nonSingWeights(h, z⃗[1 : p - 1] .- a)
        w[1 : p - 1] += c                           # Correct left weights
        w[n : -1 : n - p + 2] += c                  # Correct right weights
    end 
    w *= h                                          # Rescale weights to step size

    return w ⋅ f.(z⃗)                                # Compute integral
end
