#=
# Generalized Gregory Integration
#
# Author: Caleb Jacobs
# DLM: August 30, 2023
=#

using SpecialFunctions
using LinearAlgebra
using FFTW

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
    RHS = [-zeta(α + 1 - i) * h^i for i ∈ 0 : length(z⃗) - 1]

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
    w = getCorrections(-1, h, z⃗)   # Generate most of the weights
    w[1] -= 1                       # Fix first weight

    return w
end

"""
    nonSingInt(f, a, b, n, p)
Integrate a function `f` from `a` to `b` with no singularities using `n` nodes and order `p`
"""
function nonSingGregoryInt(f, a, b, n, p)
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

"""
    getNodes(N, r)
Generate `N` equispaced nodes on a complex circle of radius `r`.
"""
function getNodes(N::Integer, r)
    z⃗ = [r * cispi(2i / N) for i ∈ 0 : N - 1]      # Compute scaled roots of unitny

    return z⃗
end

"""
    dftζ⃗(α, h, r, N)
Construct DFT compatible zeta vector for end corrected system using `N` roots of unity scaled to have radii `r`.
"""
function dftζ⃗(α, h, r, N)
    RHS = [-zeta(α - i) * (h / r)^i for i ∈ 0 : N - 1]

    return RHS
end

"""
    getDFTCorrections(α, h, r, N)
Compute trapezoidal rule end corrections given the nodes are `N` roots of unity scaled to a radius of `r`. Note, `α` is the order of the singularity and `h` is the trapezoidal rule spacing.
"""
function getDFTCorrections(α, h, r, N)
    w⃗ = ifft(dftζ⃗(α, h, r, N))      # Compute weights using by applying an FFT to the ζ vector

    return w⃗
end

"""
    nonSingDFTInt(f, a, b, n, r, N)
Integrate a function `f` from `a` to `b` with no singularities using `n` internal nodes and `N` correction nodes a radius of `r` away from the base.
"""
function nonSingDFTInt(f, a, b, n, N, r)
    z⃗Int = [range(a, b, length = n)...]         # Internal nodes
    h = z⃗Int[2] - z⃗Int[1]                       # Internal spacing
    z⃗End = getNodes(N, r * h)                   # Correction nodes

    fInt   =  h * sum(f.(z⃗Int[2 : end - 1]))    # Internal trapezoidal rule

    if N > 0
        w = getDFTCorrections(0, h, r*h, N)     # Correction weights
        fLeft  = h * w ⋅ f.( z⃗End .+ a)         # Left end correction
        fRight = h * w ⋅ f.(-z⃗End .+ b)         # Right end correction
    else
        fLeft  = h * f(a)
        fRight = h * f(b)
    end

    ∫fdx = fLeft + fInt + fRight                # Compute corrected integral

    return ∫fdx
end

"""
    singDFTInt(f, a, b, n, r, N)
Integrate a function `f`/(x - `a`)`ᵅ` from `a` to `b` using `n` internal nodes and `N` correction nodes a radius of `r` away from the base.
"""
function singDFTInt(f, a, b, n, N, r, α)
    z⃗Int = [range(a, b, length = n)...]                         # Internal nodes
    h = z⃗Int[2] - z⃗Int[1]                                       # Internal spacing
    z⃗Ends = getNodes(N, r * h)                                  # Correction nodes

    sf(x) = f(x) / (x - a)^α                                    # Singular function

    fInt   =  h * sum(sf.(z⃗Int[2 : end - 1]))                   # Internal trapezoidal rule

    if N > 0
        wL = getDFTCorrections(α, h, r*h, N)                    # Left, singularity correction weights
        fLeft  = h^(1 - α) * wL ⋅ f.(z⃗Ends .+ a)                # Left end correction

        wR = getDFTCorrections(0, h, r*h, N)                    # Right correction weights
        fRight = h * wR ⋅ sf.(-z⃗Ends .+ b)                      # Right end correction
    else
        fLeft  = 0
        fRight = h * f(b)
    end

    ∫fdx = fLeft + fInt + fRight                                # Compute corrected integral

    return ∫fdx
end
