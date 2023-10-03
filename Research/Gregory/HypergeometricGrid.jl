#=
# Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
#
# Author: Caleb Jacobs
# DLM: October 3, 2023
=#

using SpecialFunctions
include("GenGreg.jl")

"0F0 Hypergeometric function given by ℯ^z."
function _0F0(g)
    f = exp.(g.z)

    return f
end

"1F0 Hypergeometric function given by (1 - z)^(-a)."
function _1F0(a, g)
    f = (1 .- g.z).^(-a)

    return f
end

function pFq(a, b; r = 1, n = 20, np = 3, Tr = 0.5)
    g = getGrid(n, r, ir = Tr, np = np, nl = length(b))     # Initial grid

    aIdx = length(a)                                        # Index for a
    bIdx = length(b)                                        # Index for b

    if length(a) == length(b)
        f = _0F0(g)                                         # Initial function is 0F0
    elseif length(a) == length(b) + 1
        f = _1F0(a[aIdx], g)                                # Initial funciton is 1F0
        aIdx -= 1                                           # Move a index to next layer
    else
        throw(error("Non-supported pFq order"))
    end

    for nl = length(b) : -1 : 1
        α = a[aIdx] - 1                                     # Base point singularity order
        β = b[bIdx] - a[aIdx] - 1                           # Evaluation point singularity order

        D = getDiffMat(n, r, α = α, β = β,                  # Generate differentiation matrix
                       ir = Tr, np = np, nl= nl)
        
        g = getGrid(n, r, ir = Tr, np = np, nl = nl - 1)    # Get next computation grid

        Γ = g.z.^(1 - b[bIdx]) * (                          # Front coefficient
                gamma(b[bIdx]) / 
                gamma(a[aIdx]) / 
                gamma(b[bIdx] - a[aIdx]))

        aIdx -= 1                                           # Move to next layer in a
        bIdx -= 1                                           # Move to next layer in b

        f = Γ .* D * f                                      # Compute values for next layer

        f[g.c] = 1 + 0im
    end
    
    return (g.z, f)
end
