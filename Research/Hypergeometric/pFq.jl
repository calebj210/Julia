# Compute generalized hypergeometric pFq 
#
# Author: Caleb Jacobs, Cecile Piret
# DLM: 23-02-2022

using LinearAlgebra
using SparseArrays
using SpecialFunctions

# Compute weights of center stencil
function centralWeights(α, β, z₀, zᵥ)
    N = length(zᵥ)

    Z = repeat(zᵥ, outer = (1, N))
    Z[:, 1] .= 1
    Z = cumprod(Z, dims = 2)

    ei(zₖ, pow) = zₖ^(pow + β + α + 1) * gamma(1 + β) * gamma(pow + α + 1) / gamma(pow + β + α + 2)

    return [ei(z₀[i], j) for i ∈ 1:N, j ∈ 0:(N - 1)] / Z
end

# Compute weights of endpoint correction stencil
function endStencil(α, d, p, zₖ)
    N = length(zₖ)

    order = [(0 : N - 1)...]

    A = zₖ' .^ [(0 : N - 1)...]
    vS = d .^ order .* (zeta.(-α .- order) - sum(1 ./ [(1 : p - 1)...]' .^(-α .- order), dims = 2))

    cS = A \ vS

    return cS[:]
end

# Compute pFq reduction matrix to base function
function _pFq(f, c, d, n, h, singInfo)
    NPdingNodes = 10
    singPos = Inf
    α = c[end] - 1
    β = d[end] - c[end] - 1;
    nhProb = 11                 # Problematic region node count

    idx = [(n[1] - nPdingNodes : n[2] + nPdingNodes)...]

    reIdx = idx'                # Real indices
    imIdx = idx                 # Imaginary indices
    zIdx  = reIdx .* im * imIdx # Complex compuational grid

    xᵣ = h * reIdx              # Real components of x
    xᵢ = h * imIdx              # Imaginary components of x
    z  = xᵣ + xᵢ                # Complex spatial grid points

    idxVec = zIdx[:]            # Vector of computation grid indices
    zVec   = z[:]               # Vector of spatial grid points

    F = f.(z)                   # Evaluate function over grid

    # Handle singularity info
    if isempty(signInfo) == 0 && singInfo[2] < 0
        singPos = singInfo[1]
        F[z .== singPos] .= 0
    end

    zₐ = 0                      # Base point
    IbVals = (z .- zₐ) .^ α
    IbVals[z .== zₐ] .= 0

    # Internal node indices
    @. intIdx = real(idxVec) >= n[1] && real(idxVec) <= n[2] && 
                imag(idxVec) >= n[3] && imag(idxVec) <= n[4]

    # Global operator matrices
    G1 = spzeros(length(intIdx), length(zVec)) # Standard corrections
    G2 = spzeros(length(intIdx), length(zVec)) # Standard corrections
    G  = spzeros(length(intIdx), length(zVec)) # Standard corrections
    Ga = spzeros(length(intIdx), length(zVec)) # Singularity at za = 0
    Gb = spzeros(length(intIdx), length(zVec)) # Singularity at zb = z₀
    Gc = spzeros(length(intIdx), length(zVec)) # Center corrections

    D = zeros(length(zVec))                    # Differentiation matrix

    # Endpoint corrections
    nS = 25                                    # Number of nodes in correction stencil
    t = range(0, 2π * (1 - 1 / ns), length = ns)

    Ra = 3                                     # Radius of stencil around za
    Rb = 3                                     # Radius of stencil around zb
    Rc = 3                                     # Radius of stencil around corner
    pa = 2                                     # Node away from za to start trapezoidal rule
    pa = 2                                     # Node away from zb to start trapezoidal rule
    pa = 2                                     # Node away from corner to start trapezoidal rule

    ###### Need to fix these index finders
    aIdx = findall(x -> x == round.(cis(Ra * t)), idxVec) # Grid indices of near circle elements
    bIdx = findall(x -> x == round.(cis(Rb * t)), idxVec) # Grid indices of near circle elements
    cIdx = findall(x -> x == round.(cis(Rc * t)), idxVec) # Grid indices of near circle elements
    ######

    na = length(aIdx)           # Number of nodes at orgin
    nb = length(bIdx)           # Number of nodes about evaluation point
    nc = length(cIdx)           # Number of singularity free nodes

    aIdxVec = idxVec(aIdx)
    bIdxVec = idxVec(bIdx)
    cIdxVec = idxVec(cIdx)

    # Singularity free correction stencils
    cwLft = centralWeights(0,  1,  pc, cIdxVec)
    cwRht = centralWeights(0, -1,  pc, cIdxVec)
    cwTop = centralWeights(0, -im, pc, cIdxVec)
    cwBot = centralWeights(0,  im, pc, cIdxVec)

    # Origin corrections
    cwaLft = centralWeights(α,  1,  pa, aIdxVec)
    cwaRht = centralWeights(α, -1,  pa, aIdxVec)
    cwaTop = centralWeights(α, -im, pa, aIdxVec)
    cwaBot = centralWeights(α,  im, pa, aIdxVec)

    # Evaluation point corrections
    cwbLft = centralWeights(β,  1,  pa, aIdxVec)
    cwbRht = centralWeights(β, -1,  pa, aIdxVec)
    cwbTop = centralWeights(β, -im, pa, aIdxVec)
    cwbBot = centralWeights(β,  im, pa, aIdxVec)

    r = min(nhProb, singPos / h)                                 # Radius about central node
    centerIdxVec = findall(x -> x == round.(cis(r * t)), idxVec) # Grid indices of near circle elements
    radius = r * h

    ## Stopped at line 110
end
