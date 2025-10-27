#=
#   Basic series for 2F1
#
# Author: Caleb Jacobs
# DLM: October 7, 2025
=#

function maclaurin_2f1(a, b, c, z; mord = 1000, rtol = eps())
    S = zn = coeff = 1.0

    crit_flag = false
    for n ∈ 1:mord
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        zn *= z
        val = coeff * zn
        if abs(val / S) <= rtol / 2
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

function conformal_2f1(a, b, c, z; rtol = eps(), mord = 1000)
    old = zeros(promote_type(typeof.((a,b,c))...), 4)
    old[1] = a * b / c
    for n ∈ 2:4
        old[n] = old[n-1] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end

    α = [
        1, 
         4old[1], 
        -6old[1] + 16old[2],
         4old[1] - 48old[2] +  64old[3],
    ]

    z = 1 - (1 - z)^(1/4)
    zn = z^3
    S = sum(α .* z.^(0:3))

    A(n) = -10n*(n-1)-4n*(4a+4b+1)-16a*b
    B(n) = 10(n-1)*(n-2)+6(n-1)*(4a+4b+1)+3*16a*b
    C(n) = -5(n-2)*(n-3)-4(n-2)*(4a+4b+1)-3*16a*b
    D(n) = (n-3)*(n-4)+(n-3)*(4a+4b+1)+16a*b
    E(n) = 4(n+1)*(c+n)

    crit_flag = false
    for n ∈ 3:mord
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
