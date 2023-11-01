#=
# Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
#
# Author: Caleb Jacobs
# DLM: October 31, 2023
=#

using SpecialFunctions
include("GenGreg.jl")

"Modified sign function to return 1 when z = 0"
sgn(z) = iszero(z) ? one(z) : sign(z)

"Compute roots given a power α and a branch cut rotation of θ."
function zα(z, α::Real, branch; θ::Real = 0 )
    if !branch
        if imag(z) == 0 && real(z) < 0
            return z^α * cispi(2(α % 1))
        else
            return z^α
        end
    elseif imag(z) == 0 && real(z) < 0
        return z^α * cispi(-2(α % 1))
    end

    if θ >= 0
        return angle(z) >= θ - π ? z^α : z^α * cispi( 2(α % 1))
    else
        return angle(z) <  θ + π ? z^α : z^α * cispi(-2(α % 1))
    end
end

"0F0 Hypergeometric function given by ℯᶻ."
function _0F0(g)
    f = exp.(g.z)

    return f
end

"1F0 Hypergeometric function given by (1 - z)^(-a)."
function _1F0(a, g; branch = false)
    f = (z -> zα(1 - z, -a, branch, θ = sgn.(imag(z)) * π)).(g.z)

    return f
end

"0F1 Hypergeomtreic function given by z^(1/2 - b/2) Γ(b) I_(b-1)(2√z)"
function _0F1(b, g)
    f = g.z.^((1 - b) / 2) .* gamma(b) .* besseli.(b - 1, 2sqrt.(g.z))

    return f
end

function pFq(a, b; r = 1, n = 20, np = 3, Tr = 0.5)
    o = min(length(a), length(b))                                   # Order of pFq

    g = getGrid(n, r, ir = Tr, np = np, nl = o)                     # Initial grid

    aIdx = length(a)                                                # Index for a
    bIdx = length(b)                                                # Index for b

    a = sort(a, by = abs, rev = true)                               # Sort a by modulus to help with stability
    b = sort(b, by = abs, rev = true)                               # Sort b by modulus to help with stability

    branch = false                                                  # Initialize extra branch cut correction
    if length(a) == length(b)
        f = _0F0(g)                                                 # 0F0 base function case
    elseif length(a) == length(b) + 1
        f = _1F0(a[end], g)                                         # 1F0 base function case
        h = _1F0(a[end], g, branch = true)                          # 1F0 base function case alternate branch

        branch = true

        aIdx -= 1                                                   # Move a index to next layer
    elseif length(a) == length(b) - 1
        f = _0F1(b[bIdx], g)                                        # 0F1 base function case
        f[g.c] = 1                                                  # Fix origin singularity

        bIdx -= 1                                                   # Move b index to next layer
    else
        throw(error("Non-supported pFq order"))
    end

    for nl = o : -1 : 1
        α = a[aIdx] - 1                                             # Base point singularity order
        β = b[bIdx] - a[aIdx] - 1                                   # Evaluation point singularity order

        if !branch
            D = getDiffMat(n, r, α = α, β = β,                      # Generate differentiation matrix
                       ir = Tr, np = np, nl= nl)
        else
            D0,D1,D2,D3 = getDiffMat(n, r, α = α, β = β,            # Generate differentiation matrices with alternate branch
                        ir = Tr, np = np, nl= nl, branch = branch)
        end

        g = getGrid(n, r, ir = Tr, np = np, nl = nl - 1)            # Get next computation grid

        Γ = g.z.^(1 - b[bIdx]) * (                                  # Front coefficient
                gamma(b[bIdx]) / 
                gamma(a[aIdx]) / 
                gamma(b[bIdx] - a[aIdx]))

        aIdx -= 1                                                   # Move to next layer in a
        bIdx -= 1                                                   # Move to next layer in b

        if !branch
            f = Γ .* D * f                                          # Compute values for next layer
        else
            f, h = (Γ .* (D0 * f + D1 * h), Γ .* (D2 * f + D3 * h)) # Compute values for next layer
        end

        f[g.c] = 1 + 0im
    end
    
    return (g.z, f)
end
