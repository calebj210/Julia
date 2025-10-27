#=
#   Routine for computing 2F1 by comparing multiple evaluations
#
# Author: Caleb Jacobs
# DLM: October 27, 2025
=#

include("Transformations.jl")

function comparison_2f1(a, b, c, z; rtol = 1e-13, ord = 4)
    if real(z) <= 0.5
        if abs2(z) < 1
            trans = [:z, :zoverzminusone]
        else
            trans = [:zoverzminusone, :oneoverz, :oneoveroneminusz]
        end
    else
        if abs2(1 - z) < 1
            trans = [:oneminusz, :oneminusoneoverz]
        else
            trans = [:oneminusoneoverz, :oneoverz, :oneoveroneminusz]
        end
    end

    val = compare(a, b, c, z, trans; rtol = rtol, ord = ord)
    return val
end

function compare(a, b, c, z, trans; ord = 4, kwargs...)
    dif = Inf
    idx = [0,0]

    N = length(trans)
    vals = zeros(ComplexF64, ord*N)
    for n ∈ 1:ord*N - 1
        ordn = ord - ((n - 1) ÷ N)
        if iszero(vals[n])
            vals[n] = transformations[((n - 1) % N) + 1](a, b, c, z, ordn)
        end
        for k ∈ (n+1):ord*N
            ordk = ord - ((k - 1) ÷ N)
            if iszero(vals[k])
                vals[k] = transformations[((k - 1) % N) + 1](a, b, c, z, ordk)
            end

            if isapprox(vals[n], vals[k]; kwargs...)
                return (vals[n] + vals[k]) / 2
            elseif abs(vals[n] - vals[k]) < dif
                dif = abs(vals[n] - vals[k])
                idx = [n,k]
            end
        end
    end

    @warn "Tolerance not met, answer within a relative tolerance of $(abs(dif / vals[first(idx)]))."

    return vals[first(idx)]
end
