#=
# Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
#
# Author: Caleb Jacobs
# DLM: October 31, 2024
=#

using SpecialFunctions
include("GenGreg.jl")
include("Z1Expansion.jl")

"Modified sign function to return 1 when z = 0"
sgn(z)  = iszero(z) ?  one(z) : sign(z)
nsgn(z) = iszero(z) ? -one(z) : sign(z)

"Compute roots given a power α and a branch cut rotation of θ."
function oneMinusZα(z, α::Real, branch = false)
    if !branch
        if imag(z) == 0 && real(z) > 1
            return (1 - z)^α * cispi(2(α % 1))
        else
            return (1 - z)^α
        end
    else
        if imag(z) == 0 && real(z) > 1
            return (1 - z)^α
        elseif imag(z) > 0
            return (1 - z)^α * cispi(2(α % 1))
        else
            return (1 - z)^α * cispi(-2(α % 1))
        end
    end
end

"0F0 Hypergeometric function given by ℯᶻ."
function _0F0(g)
    f = exp.(g.z)

    return f
end

"1F0 Hypergeometric function given by (1 - z)^(-a)."
function _1F0(a, g; branch = false)
    f = oneMinusZα.(g.z, -a, branch)

    return f
end

"0F1 Hypergeomtreic function given by z^(1/2 - b/2) Γ(b) I_(b-1)(2√z)"
function _0F1(b, g)
    f = g.z.^((1 - b) / 2) .* gamma(b) .* besseli.(b - 1, 2sqrt.(g.z))

    return f
end

function pFq(a, b; grid_radius = 1.99, grid_points = 41, padding_layers = 5, taylor_radius = 0.5, modify_z1 = true, correction_radius = .5, inner_radius = .6, outer_radius = .8, z1_expansion_order = 70)
    o = min(length(a), length(b))                                   # Order of pFq

    g = getGrid(grid_points, grid_radius, ir = taylor_radius, np = padding_layers, nl = o)                     # Initial grid

    aIdx = length(a)                                                # Index for a
    bIdx = length(b)                                                # Index for b

#     a = sort(a, by = real, rev = true)                              # Sort a by modulus to help with stability
#     b = sort(b, by = real, rev = true)                              # Sort b by modulus to help with stability

    branch = false                                                  # Initialize extra branch cut correction
    if length(a) == length(b)
        f = _0F0(g)                                                 # 0F0 base function case
        h = f
    elseif length(a) == length(b) + 1
        f = _1F0(a[end], g)                                         # 1F0 base function case
        h = _1F0(a[end], g, branch = true)                          # 1F0 base function case alternate branch
        ωa = [1]                                                    # Initial z = 1 expansion coefficient

        branch = true

        aIdx -= 1                                                   # Move a index to next layer
    elseif length(a) == length(b) - 1
        f = _0F1(b[bIdx], g)                                        # 0F1 base function case
        f[g.c] = 1                                                  # Fix origin singularity
        h = f

        bIdx -= 1                                                   # Move b index to next layer
    else
        throw(error("Non-supported pFq order"))
    end

    for nl = o : -1 : 1
        α = a[aIdx] - 1                                             # Base point singularity order
        β = b[bIdx] - a[aIdx] - 1                                   # Evaluation point singularity order

        if !branch
            D0, D1 = getDiffMat(grid_points, grid_radius, α = α, β = β,                 # Generate differentiation matrix
                       ir = taylor_radius, np = padding_layers, nl= nl)
        else
            D0,D1,D2,D3 = getDiffMat(grid_points, grid_radius, α = α, β = β,            # Generate differentiation matrices with alternate branch
                        ir = taylor_radius, np = padding_layers, nl= nl, branch = branch)
        end

        g = getGrid(grid_points, grid_radius, ir = taylor_radius, np = padding_layers, nl = nl - 1)            # Get next computation grid

        Γ = g.z.^(1 - b[bIdx]) * (                                  # Front coefficient
                gamma(b[bIdx]) / 
                gamma(a[aIdx]) / 
                gamma(b[bIdx] - a[aIdx]))

        if !branch
            f = Γ .* (D0 * f + D1 * h)                              # Compute values for next layer
            h = f
        else
            f, h = (Γ .* (D0 * f + D1 * h), 
                    Γ .* (D2 * f + D3 * h))                         # Compute values for next layer

            if modify_z1
                zm1 = abs.(g.z .- 1)
                cIdx = zm1 .<= correction_radius                                # Corrected node indices

                zIdx = (inner_radius .<= zm1) .&& (zm1 .< outer_radius)         # Stencil nodes
                zCirc = g.z[zIdx]
                fCirc = f[zIdx]

                (ωa, ωb) = getZ1ExpansionWeights(a[aIdx : end], 
                                                 b[bIdx : end], 
                                                 ωa, zCirc, fCirc,
                                                 n = z1_expansion_order)

                f[cIdx] = z1PFQ(a[aIdx : end], b[bIdx : end], ωa, ωb, g.z[cIdx])
                h[cIdx] = z1PFQ(a[aIdx : end], b[bIdx : end], ωa, ωb, g.z[cIdx], branch = true)
            end
        end

        aIdx -= 1                                                   # Move to next layer in a
        bIdx -= 1                                                   # Move to next layer in b
        f[g.c] = 1 + 0im                                            # Fix origin
    end
    
    if !branch
        return (g.z, f)
    else
        return (g.z, f, h)
    end
end
