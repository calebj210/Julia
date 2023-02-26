# Compute generalized hypergeometric pFq 
#
# Authors: Caleb Jacobs, Cecile Piret
# DLM: 26-02-2022

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
    @. intIdx = findall(real(idxVec) >= n[1] && real(idxVec) <= n[2] && 
                        imag(idxVec) >= n[3] && imag(idxVec) <= n[4])

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

    Ra = 3                                      # Radius of stencil around za
    Rb = 3                                      # Radius of stencil around zb
    Rc = 3                                      # Radius of stencil around corner
    pa = 2                                      # Node away from za to start trapezoidal rule
    pa = 2                                      # Node away from zb to start trapezoidal rule
    pa = 2                                      # Node away from corner to start trapezoidal rule

    # Find grid points that are approximately on circle
    aIdx = findall(x -> !isnothing(x), 
           indexin(idxVec, round.(Ra * cis.(t))))
    bIdx = findall(x -> !isnothing(x), 
           indexin(idxVec, round.(Rb * cis.(t))))
    cIdx = findall(x -> !isnothing(x), 
           indexin(idxVec, round.(Rc * cis.(t))))

    na = length(aIdx)                           # Number of nodes at orgin
    nb = length(bIdx)                           # Number of nodes about evaluation point
    nc = length(cIdx)                           # Number of singularity free nodes

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

    r = min(nhProb, singPos / h)                    # Radius about central node
    originIdx = findall(x -> !isnothing(x),         # origin index
           indexin(idxVec, round.(r * cis.(t))))
    radius = r * h

    origin = findall(iszero, zIdx)                  # Origin position
    
    for i ∈ eachindex(intIdx)
        Wtr = zeros(size(z))                        # Trapezoidal rule weights
        Wa  = zeros(length(z))
        Wb  = zeros(length(z))
        Wc  = zeros(length(z))

        zb = zVec(intIdx(i))                        # Evaluation point

        if abs(zb) <= radius
            Gc[i, originIdx] = centralWeights(α, β, zb, idxVec[originIdx]) 
        else
            x0Idx = round(real(zb) / h)             # real part of evaluation point
            y0Idx = round(imag(zb) / h)             # imaginary part of evaluation point
            
            δza = 10                                # Safe distance from za
            δzb = 10                                # Safe distance from zb
            δzc = 0                                 # Safe distance from zc
            
            flagb = 0                               # Relocate branch cut
            flagn = 0                               # Allow integration along negative real axis
            if x0Idx == 0 || y0Idx == 0
                flagb = 0
                path = [0 0; x0Idx y0Idx]
            elseif x0Idx >= -δza && abs(y0Idx) > δzb
                flagb = 1
                if y0Idx < 0
                    falgn = 1
                end
                path = [0 0; -(δza + δzb) 0; -(δza + δzb) y0Idx; x0Idx y0Idx]
            elseif x0Idx >= 0 && abs(y0Idx) < δzb
                falgb = 1
                path = [0 0; 0 sign(y0Idx)*(δza + δzb); x0Idx sign(y0Idx)*(δza + δzb); x0Idx y0Idx]
            elseif x0Idx < 0 && y0Idx > δzb         # Positive test
                flagb = 0
                path = [0 0; x0Idx 0; x0Idx y0Idx]
            elseif x0Idx < 0 && y0Idx < -δzb        # Negative test
                flagn = 1
                flagb = 0
                path = [0 0; x0Idx 0; x0Idx y0Idx]
            elseif x0Idx < 0 && abs(y0Idx) <= δzb
                flagb = 0
                path = [0 0; 0 sign(y0Idx)*(δza + δzb); x0Idx sign(y0Idx)*(δza + δzb); x0Idx y0Idx]
            end

            IaVals = similar(z)
            if flagb == 1
                IaVals .= (zb - z).^β
            else
                IaVals .= (z - zb).^β
            end
            IabVals = IaVals .* IbVals

            # Integrate over path
            for pathIdx ∈ 2 : size(path, 1)
                # Get TR weights
                stpx = sign(path[pathIdx, 1] - path[pathIdx - 1, 1])
                dx   = stpx + (stpx == 0)

                stpy = sign(path[pathIdx, 2] - path[pathIdx - 1, 2])
                dy   = stpy + (spty == 0)

                yiPath = origin[1][1] .- [path[pathIdx - 1, 2] : dy : path[pathIdx, 2]...]
                xiPath = origin[1][2] .- [path[pathIdx - 1, 1] : dx : path[pathIdx, 1]...]

                dir = stpx + im * stpy
                dir = dir / abs(dir)
                hDir = h * dir

                corrTR = ones(size(z[yiPath, xiPath]))

                if flagn == 1 && abs(zIdx[yiPath[1], xiPath[1]]) < eps() && abs(imag(zIdx[yiPath[end], xiPath[end]])) < eps()
                    corrTR = corrTR * cis(-2π * α)
                end

                trPath = hDir * corrTR

                if pathIdx == 2 && pathIdx == size(path, 1)
                    trPath[1 : pa] .= 0
                    trPath[end + 1 - pb : end] .= 0
                elseif pathIdx == 2
                    trPath[1 : pa] .= 0
                    trPath[end + 1 - pc : end] .= 0
                elseif pathIdx == size(path, 1)
                    trPath[1 : pc] .= 0
                    trPath[end + 1 - pb : end] .= 0
                else 
                    trPath[1 : pc] .= 0
                    trPath[end + 1 - pc : end] .= 0
                end

                Wtr[yiPath, xiPath] += trPath       # TR weights

                # Get start corrections
                if pathIdx == 2
                    zsaImag = h * imag(aIdxVec)

                    sheetCorra = ones(na, 1)
                    if imag(zb) <= 2 * Ra * h > 0 && real(zb) >= 0 && flagb == 0
                        sheetCorra[zsaImag .>= imag(zb)] .= cis(-2π * β)
                    elseif imag(zb) >= -2 * Ra * h && imag(zb) <= 0 && real(zb) >= 0 && flagb == 0
                        sheetCorra[zsaImag .< imag(zb)] .= cis(2π * β)
                    end

                    if flagn == 1
                        sheetCorra *= cis(-2π * α)
                    end

                    Wa[aIdx] += hDir^(α + 1) * sheetCorra .*(
                                (spty == 0 && stpx >  0) * cwaLft + # Left start correction
                                (stpy == 0 && stpx <  0) * cwaRht + # Left end correction
                                (stpy >  0 && stpx == 0) * cwaBot + # Bottom start correction
                                (stpy <  0 && stpx == 0) * cwaTop)  # Top start correction
                else
                    cjIdx = findall(indexin(idxVec, zIdx[yiPath[1], xiPath[1]] + cIdxVec))
                    Zc    = zIdx[yiPath[1], xiPath[1]]

                    sheetCorrc = ones(nc, 1)
                    if imag(Zc) <= 2Rc && imag(Zc) >= 0 && real(Zc) < 0 && flagn == 0
                        sheetCorrc[imag[idxVec[cjIdx]] .< 0] .= cis(2π * α)
                    elseif (imag(Zc) >= -2Rb && imag(Zc) < 0 && real(Zc) < 0 && flagn == 0) || (imag(Zc) == 0 && real(Zc) < 0 && flagn == 1)
                        sheetCorrc[imag[idxVec[cjIdx]] .< 0] .= cis(-2π * α)
                    end

                    Wc[cIdx] += hDir * sheetCorrc .* (
                                (spty == 0 && stpx >  0) * cwLft +  # Left start correction
                                (stpy == 0 && stpx <  0) * cwRht +  # Left end correction
                                (stpy >  0 && stpx == 0) * cwBot +  # Bottom start correction
                                (stpy <  0 && stpx == 0) * cwTop)   # Top start correction
                end

                # End correction
                if pathIdx == size(path, 1)
                    trnsltdbIdxVec = zIdx[yiPath[end], xiPath[end]] + bIdxVec
                    bIdx = findall(indexin(idxVec, trnsltdbIdxVec))
                    zsbImag = imag(trnsltdbIdxVec)

                    sheetCorrb = ones(nb, 1)
                    if imag(zb) <= Rb * h && imag(zb) >= 0 && real(zb) < 0
                        sheetCorrb[zsbImag .< 0] .= cis(2π * α)
                    elseif imag(zb) >= -Rb * h && imag(zb) < 0 && real(zb) < 0
                        sheetCorrb[zsbImag .>= 0] .= cis(-2π * α)
                    elseif imag(zb) >= -2Rb * h && imag(zb) < 0 && real(zb) > 1 && length(c) == 2
                        sheetCorrb[zsbImag .> 0] .= cis(2π * c[1])
                    elseif imag(zb) <= 2Rb * h && imag(zb) > 0 && real(zb) > 1 && length(c) == 2
                        sheetCorrb[zsbImag .<= 0] .= cis(-2π * c[1])
                    end

                    Wb[bIdx] += (-hDir) ^ (β + 1) * sheetCorrb .* (
                                (spty == 0 && stpx >  0) * cwbLft +  # Right end correction
                                (stpy == 0 && stpx <  0) * cwbRht +  # Left end correction
                                (stpy >  0 && stpx == 0) * cwbBot +  # Bottom end correction
                                (stpy <  0 && stpx == 0) * cwbTop)   # Top end correction

                else
                    cjIdx = findall(indexin(idxVec, zIdx[yiPath[end], xiPath[end]] + cIdxVec))
                    Zc    = zIdx[yiPath[end], xiPath[end]]

                    sheetCorrc = ones(nc, 1)
                    if imag(Zc) <= 2Rc && imag(Zc) >= 0 && real(Zc) < 0 && flagn == 0
                        sheetCorrc[imag[idxVec[cjIdx]] .< 0] .= cis(2π * α)
                    elseif (imag(Zc) >= -2Rb && imag(Zc) < 0 && real(Zc) < 0) || (imag(Zc) == 0 && real(Zc) < 0 && flagn == 1)
                        sheetCorrc[imag[idxVec[cjIdx]] .< 0] .= cis(-2π * α)
                    end

                    Wc[cIdx] += hDir * sheetCorrc .* (
                                (spty == 0 && stpx >  0) * cwLft +  # Left start correction
                                (stpy == 0 && stpx <  0) * cwRht +  # Left end correction
                                (stpy <  0 && stpx == 0) * cwBot +  # Bottom start correction
                                (stpy >  0 && stpx == 0) * cwTop)   # Top start correction
                end
            end

            # Transform the matrices W with Wb into the jth row of the global matrix G
            tmp = (Wc .* IabVals[:]).'
            Gpos1 = findall(tmp)
            G1[intIdx[i], Gpos1] .= tmp[Gpos1]

            tmp = (Wtr[:] .* IabVals[:]).'
            Gpos2 = findall(tmp)
            G2[intIdx[i], Gpos2] .= tmp[Gpos2]

            tmp = ((Wc + Wtr[:]) .* IabVals[:]).'
            Gpos = findall(tmp)
            G[intIdx[i], Gpos] .= tmp[Gpos]

            tmp = (Wb .* IbVals[:]).'
            Gposb = findall(tmp)
            G1[intIdx[i], Gposb] .= tmp[Gposb]

            tmp = (Wa .* IabVals[:]).'
            Gposa = findall(tmp)
            G2[intIdx[i], Gposa] .= tmp[Gposa]
        end

        if angle(zb) > 0 && flagb == 0
            Cb[i] = cis(π * β)
            Ca[i] = cis(π * β)
        elseif flagb == 0
            Cb[i] = cis(-π * β)
            Ca[i] = cis(-π * β)
        end
        
        if flagb == 1
            if y0Idx < 0 && y0Idx >= -10
                Cb[i] = cis(π * β)
            elseif y0Idx <= -10
                Cb[i] = cis(-π * β)
            elseif y0Idx > 0 && y0Idx <= 10
                Cb[i] = cis(-sign(imag(zb)) * π * β)
            elseif y0Idx > 10
                Cb[i] = cis(-sign(imag(zb)) * π * β)
            end

            Ca[i] = 1
        end
    end

    Dmat = GCenter + diag(Cb) * Gb[intIdx, :] - diag(Ca) * (G1[intIdx, :] + Ga[intIdx, :])
    int = Dmat * F[:]

    return gamma(d[end]) / gamma(c[end]) / gamma(d[end] - c[end]) * int .* z[] .^ (1 - d[end])
end
