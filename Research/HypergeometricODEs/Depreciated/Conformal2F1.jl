#=
# Compute 2F1 via transformations and a conformal mapping
#
# Author: Caleb Jacobs
# DLM: September 26, 2025
=#

using DSP
using LinearAlgebra
using SpecialFunctions

include("Transformations.jl")

# Generate series weights for z -> 1 - (1 - ζ)^deg
function conformalweights(ord; deg = 4)
    if deg == 2
        a = [0, 2, -2]
    elseif deg == 3
        a = [0, 3, -3, 1]
    elseif deg == 4
        a = [0, 4, -6, 4, -1]
    elseif deg == 5
        a = [0, 5, -10, 10, -5, 1]
    elseif deg == 6
        a = [0, 6, -15, 20, -15, 6, -1]
    else
        error("Degree deg = $(deg) not supported.")
    end
    A = zeros(Float64, ord, ord)
    A[1,1] = 1
    for i ∈ 2:ord
        A[:,i] = conv(a,A[:,i - 1])[1:ord]
    end

    return LowerTriangular(A)
end

function conformal_2f1(a, b, c, z; rtol = eps(), mord = 1000, raw = false)
    old = zeros(4)
    old[1] = a * b / c
    for n ∈ 2:4
        old[n] = old[n-1] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end

    α = [
        1, 
         4old[1], 
        -6old[1] + 16old[2],
         4old[1] - 48old[2] +  64old[3],
         -old[1] + 68old[2] - 288old[3] + 256old[4]
    ]

    # α = conformalweights(5) * old

    if !raw
        z = 1 - (1 - z)^.25
    end

    if abs2(z) >= 1
        return [NaN + im * NaN, NaN + im * NaN]
    end

    zn = z^4
    S = dot(conj(z .^ (0:4)), α)

    A(n) = -10n*(n-1)-4n*(4a+4b+1)-16a*b
    B(n) = 10(n-1)*(n-2)+6(n-1)*(4a+4b+1)+3*16a*b
    C(n) = -5(n-2)*(n-3)-4(n-2)*(4a+4b+1)-3*16a*b
    D(n) = (n-3)*(n-4)+(n-3)*(4a+4b+1)+16a*b
    E(n) = 4(n+1)*(c+n)

    crit_flag = false
    for n ∈ 4:mord
        push!(α, -(A(n) * α[n + 1] + B(n) * α[n] + C(n) * α[n - 1] + D(n) * α[n - 2]) / E(n))

        zn *= z
        val = last(α) * zn
        if abs(val / S) <= rtol
            if crit_flag
                break
            else
                crit_flag = true
            end
        elseif crit_flag
            crit_flag = false
        end
        S += val
    end

    return [S, NaN + NaN * im]
end

function 
end
