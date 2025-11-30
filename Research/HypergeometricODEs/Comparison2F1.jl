#=
#   Routine for computing 2F1 by comparing multiple evaluations
#
# Author: Caleb Jacobs
# DLM: November 30, 2025
=#

include("Transformations.jl")

function comparison_2f1(a, b, c, z; rtol = 1e-14, ord = 2, esterr = false)
    if isreal(z)
        z = real(z) - 0im
    end

    if abs2(z) < 1 && abs2(1 - z) < 1
        if real(z) <= 0.5
            trans = [:z, :zoverzminusone]
        else
            trans = [:oneminusz, :oneminusoneoverz]
        end
    else
        if real(z) <= 0.5
            trans = [:oneoverz, :oneoveroneminusz, :zoverzminusone, :z, :oneminusz]
        else
            trans = [:oneoverz, :oneoveroneminusz, :oneminusoneoverz, :z, :oneminusz]
        end
    end

    # if abs2(z) < 1
    #     trans = [:z, :zoverzminusone, :oneminusz, :oneminusoneoverz]
    # else
        # trans = [:z, :oneminusz, :zoverzminusone, :oneminusoneoverz, :oneoverz, :oneoveroneminusz]
    # end

    val = compare(a, b, c, z, trans; rtol, ord, esterr)
    return val
end

function compare(a, b, c, z, trans; ord = 4, esterr = false, kwargs...)
    dif = Inf
    idx = [0,0]

    N = length(trans)
    vals = zeros(ComplexF64, ord*N)
    for n ∈ 1:ord*N - 1
        ordn = ((n - 1) ÷ N) + 1
        trann = trans[((n - 1) % N) + 1]
        if iszero(vals[n])
            vals[n] = transformations[trann](a, b, c, z, ordn)
        end
        for k ∈ (n+1):ord*N
            ordk = ((k - 1) ÷ N) + 1
            trank = trans[((k - 1) % N) + 1]
            if iszero(vals[k])
                vals[k] = transformations[trank](a, b, c, z, ordk)
            end

            if isapprox(vals[n], vals[k]; kwargs...)
                if esterr
                    return ((vals[n] + vals[k]) / 2, abs(vals[n] - vals[k]) / max(abs(vals[n]), abs(vals[k])))
                else
                    return (vals[n] + vals[k]) / 2
                end
            elseif abs(vals[n] - vals[k]) < dif
                dif = abs(vals[n] - vals[k])
                idx = [n,k]
            end
        end
    end

    # @warn "Tolerance not met, answer within a relative tolerance of $(abs(dif / vals[first(idx)]))."

    if esterr
        return (sum(vals[idx]) / 2, dif / max(abs.(vals[idx])...))
    else
        return sum(vals[idx]) / 2
    end
end
