#=
# Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
#
# Author: Caleb Jacobs
# DLM: May 21, 2024
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

function pFq(a, b; r = 1.99, n = 41, np = 5, Tr = 0.5, modifyZ1 = true, corrR = .25, circR = .8, branchN = 150, circN = 150, interpN = 10)
    o = min(length(a), length(b))                                   # Order of pFq

    g = getGrid(n, r, ir = Tr, np = np, nl = o)                     # Initial grid

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
        ωa = [1]                                                    # Initial z = 1 expansion coeficient

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
            D0, D1 = getDiffMat(n, r, α = α, β = β,                 # Generate differentiation matrix
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

        if !branch
            f = Γ .* (D0 * f + D1 * h)                              # Compute values for next layer
            h = f
        else
            f, h = (Γ .* (D0 * f + D1 * h), 
                    Γ .* (D2 * f + D3 * h))                         # Compute values for next layer

            if modifyZ1
                cIdx = abs.(g.z .- 1) .<= corrR                     # Corrected node indices

                zCirc = getZVals(r = circR, n = circN)              # Scaled roots of unity
                
                fCirc = Vector{ComplexF64}(undef, 0)                # Interpolated function values
                ω = barycentricWeights(getInterpolantNodes(g, g.z[g.c], n = interpN))
                for zc ∈ zCirc
                    zIdx = getInterpolantNodes(g, zc, n = interpN)
                    if real(zc) >= 1.0
                        fTmp = [nsgn(imag(zc)) == nsgn(imag(g.z[zi])) ? f[zi] : h[zi] for zi ∈ zIdx]
                    else
                        fTmp = f[zIdx]
                    end
                    push!(fCirc, barycentricInterpolate(g.z[zIdx], ω, fTmp)(zc))
                end

                ωaTmp = ωa
                (ωa, ωb) = getZ1ExpansionWeights(a[aIdx : end], 
                                                 b[bIdx : end], 
                                                 ωa, zCirc, fCirc,
                                                 n = branchN)

                f[cIdx] = z1PFQ(a[aIdx : end], b[bIdx : end], ωa, ωb, g.z[cIdx])
                h[cIdx] = f[cIdx] + Φ(a[aIdx : end], b[bIdx : end], ωaTmp, g.z[cIdx], n = branchN)
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
