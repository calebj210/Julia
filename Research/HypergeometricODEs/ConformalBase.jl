#=
#   Conformal mapping base functions
#
# Author: Caleb Jacobs
# DLM: October 26, 2025
=#

function recurrence_relations(a, b, c, ord)
    B = 1 + ord * (a + b)
    C = ord^2 * a * b
    if ord == 1
        return n -> [
            -(a + n) * (b + n),
            (c + n) * (n + 1)
        ]
    elseif ord == 2
        return n -> [
            C + B*(n-1) + (n-1)*(n-2),
            -C - 2B*n - 3(n+0)*(n-1),
            2c*(n+1) + 2(n+1)*(n+0)
        ]
    elseif ord == 3
        return n -> [
            -C - B*(n-2) - (n-2)*(n-3),
            2C + 3B*(n-1) + 4(n-1)*(n-2),
            -C - 3B*(n) - 6(n+0)*(n-1),
            3c*(n+1)  + 3(n+1)*n
        ]
    elseif ord == 4
        return n -> [
            (n-3)*(n-4)+(n-3)*B+C,
            -5(n-2)*(n-3)-4(n-2)*B-3*C,
            10(n-1)*(n-2)+6(n-1)*B+3*C,
            -10n*(n-1)-4n*B-C,
            4(n+1)*(c+n)
        ]
    elseif ord == 5
        return n -> [
            -  (n-4)*(n-5) -   B*(n-4) -  C,
              6(n-3)*(n-4) +  5B*(n-3) + 4C,
            -15(n-2)*(n-3) - 10B*(n-2) - 6C,
             20(n-1)*(n-2) + 10B*(n-1) + 4C,
            -15(n+0)*(n-1) -  5B*(n-0) -  C,
              5(n+1)*(n-0) +  5c*(n+1)
        ]
    elseif ord == 6
        return n -> [
            -  (n-6)*(n-7) -   B*(n-6) -   C,
              8(n-5)*(n-6) +  7B*(n-5) +  6C,
            -28(n-4)*(n-5) - 21B*(n-4) - 15C,
             56(n-3)*(n-4) + 35B*(n-3) + 20C,
            -70(n-2)*(n-3) - 35B*(n-2) - 15C,
             56(n-1)*(n-2) + 21B*(n-1) +  6C,
            -27(n+0)*(n-1) - 6(B+c)*n  -   C,
              6(n+1)*(n+0) + 6c* (n+1)
        ]
    else
        error("Unsupported conformal mapping order $(ord)")
    end
end

function initialize_series(a, b, c, ord)
    α = [one(promote_type(typeof.((a, b, c))...))]
    S = zn = one(ComplexF64)
    
    A = recurrence_relations(a, b, c, ord)

    return (α, zn, S, A)
end

function next_coefficient!(α, A)
    n = length(α) - 1
    An = A(n)
    N = length(An) - 1

    if n < N
        push!(α, -sum(An[end - (n + 1):end - 1] .* α) / last(An))
    else
        push!(α, -sum(An[1:end - 1] .* α[end - N + 1:end]) / last(An))
    end
end

function conformal_2f1(a, b, c, z, ord::Integer = 4; maxord = 100, rtol = eps())
    z = 1 - (1 - z)^(1 / ord)

    α, zn, S, A = initialize_series(a, b, c, ord)

    crit_flag = false
    for _ ∈ length(α):maxord
        next_coefficient!(α, A)

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
