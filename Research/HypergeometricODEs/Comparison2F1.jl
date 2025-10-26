#=
#   Routine for computing 2F1 by comparing multiple evaluations
#
# Author: Caleb Jacobs
# DLM: October 26, 2025
=#

include("Transformations.jl")

function comparison_2f1(a, b, c, z; rtol = 1e-13)
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

    val = compare(a, b, c, z, trans; rtol = rtol)
    return val
end

function compare(a, b, c, z, trans; kwargs...)
    N = length(trans)
    vals = zeros(ComplexF64, 4N)
    for n ∈ 1:4N - 1
        ordn = 4 - ((n - 1) ÷ N)
        if iszero(vals[n])
            vals[n] = transformations[((n - 1) % N) + 1](a, b, c, z, ordn)
        end
        for k ∈ (n+1):4N
            ordk = 4 - ((k - 1) ÷ N)
            if iszero(vals[k])
                vals[k] = transformations[((k - 1) % N) + 1](a, b, c, z, ordk)
            end

            if isapprox(vals[n], vals[k]; kwargs...)
                return vals[n]
            end
        end
    end

    return NaN * (0 + 0im)
end
