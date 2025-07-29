#=
# Compute 2F1 via transformations and a conformal mapping
#
# Author: Caleb Jacobs
# DLM: July 29, 2025
=#

using DSP
using LinearAlgebra

# Generate series weights for z -> 1 - (1 - ζ)^4
function conformalweights(ord)
    a = [0, 4, -6, 4, -1]
    A = zeros(Float64, ord, ord)
    A[1,1] = 1
    for i ∈ 2:ord
        A[:,i] = conv(a,A[:,i - 1])[1:ord]
    end

    return LowerTriangular(A)
end

function conformal_maclaurin_2f1(a, b, c, z; rtol = eps(), mord = 1000)
    old = zeros(5)
    old[1] = 1
    for n ∈ 1:4
        old[n+1] = old[n] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end

    α = conformalweights(5) * old

    z = 1 - (1 - z)^.25

    if abs2(z) >= 1
        return NaN + im * NaN
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

    return S
end

function convolutional_maclaurin_2f1(a, b, c, z; max_ord = 25, rtol = 1e-15)
    z = 1 - (1 - z)^.25
    S = zn = coeff = 1.0
    coeffs = [1.0 + 0im]

    # crit_flag = false
    for n ∈ 1:max_ord
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        # zn *= z
        # val = coeff * zn
        # if abs(val / S) <= rtol
        #     if crit_flag
        #         break
        #     else
        #         crit_flag = true
        #     end
        # elseif crit_flag
        #     crit_flag = false
        # end
        # S += val

        push!(coeffs, coeff)
    end

    return dot(conj(z .^ (0:length(coeffs) - 1)), conformalweights(length(coeffs)), coeffs)
end
