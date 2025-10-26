function convolutional_maclaurin_2f1(a, b, c, z; max_ord = 25, deg = 4)
    z = 1 - (1 - z)^(1 / deg)
    coeff = 1.0
    coeffs = [1.0 + 0im]

    for n ∈ 1:max_ord
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        push!(coeffs, coeff)
    end

    return dot(conj(z .^ (0:length(coeffs) - 1)), conformalweights(length(coeffs); deg = deg), coeffs)
end

function conformal_2_2f1(a, b, c, z; rtol = eps(), mord = 100, raw = false)
    old = zeros(2)
    old[1] = 1
    old[2] = a*b/c
    α = conformalweights(2; deg = 2) * old

    if !raw
        z = 1 - sqrt(1 - z)
    end

    if abs2(z) >= 1
        return NaN + im * NaN
    end

    zn = z
    S = dot(conj(z .^ (0:1)), α)

    B = 1 + 2(a + b)
    C = 4a*b 
    A(n) = [
        C + B*(n-1) + (n-1)*(n-2),
        -C - 2B*n - 3(n+0)*(n-1),
        2c*(n+1) + 2(n+1)*(n+0)
    ]

    crit_flag = false
    for n ∈ 1:mord
        tmp = A(n)
        push!(α, -sum(tmp[1:2] .* α[end-1:end]) / last(tmp))

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

function conformal_3_2f1(a, b, c, z; rtol = eps(), mord = 100, raw = false)
    old = zeros(3)
    old[1] = 1
    for n ∈ 1:2
        old[n+1] = old[n] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end
    α = conformalweights(3; deg = 3) * old

    if !raw
        z = 1 - (1 - z)^(1/3)
    end

    if abs2(z) >= 1
        return NaN + im * NaN
    end

    zn = z^2
    S = dot(conj(z .^ (0:2)), α)

    B = 1 + 3(a + b)
    C = 9a*b 
    A(n) = [
        -C - B*(n-2) - (n-2)*(n-3),
        2C + 3B*(n-1) + 4(n-1)*(n-2),
        -C - 3B*(n) - 6(n+0)*(n-1),
        3c*(n+1)  + 3(n+1)*n
    ]

    crit_flag = false
    for n ∈ 2:mord
        tmp = A(n)
        push!(α, -sum(tmp[1:3] .* α[end-2:end]) / last(tmp))

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

function conformal_5_2f1(a, b, c, z; rtol = eps(), mord = 100, raw = false)
    old = zeros(5)
    old[1] = 1
    for n ∈ 1:4
        old[n+1] = old[n] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end
    α = conformalweights(5; deg = 5) * old

    if !raw
        z = 1 - (1 - z)^(1/5)
    end

    if abs2(z) >= 1
        return NaN + im * NaN
    end

    zn = z^4
    S = dot(conj(z .^ (0:4)), α)

    B = 1 + 5(a + b)
    C = 25a*b
    A(n) = [
        -  (n-4)*(n-5) -   B*(n-4) -  C,
          6(n-3)*(n-4) +  5B*(n-3) + 4C,
        -15(n-2)*(n-3) - 10B*(n-2) - 6C,
         20(n-1)*(n-2) + 10B*(n-1) + 4C,
        -15(n+0)*(n-1) -  5B*(n-0) -  C,
          5(n+1)*(n-0) +  5c*(n+1)
    ]

    crit_flag = false
    for n ∈ 4:mord
        tmp = A(n)
        push!(α, -sum(tmp[1:5] .* α[end-4:end]) / last(tmp))

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

function conformal_6_2f1(a, b, c, z; rtol = eps(), mord = 100, raw = false)
    old = zeros(7)
    old[1] = 1
    for n ∈ 1:6
        old[n+1] = old[n] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end
    α = conformalweights(7; deg = 6) * old

    if !raw
        z = 1 - (1 - z)^(1/6)
    end

    if abs2(z) >= 1
        return NaN + im * NaN
    end

    zn = z^6
    S = dot(conj(z .^ (0:6)), α)

    B = 1 + 6(a + b)
    C = 36a*b
    A(n) = [
        -  (n-6)*(n-7) -   B*(n-6) -   C,
          8(n-5)*(n-6) +  7B*(n-5) +  6C,
        -28(n-4)*(n-5) - 21B*(n-4) - 15C,
         56(n-3)*(n-4) + 35B*(n-3) + 20C,
        -70(n-2)*(n-3) - 35B*(n-2) - 15C,
         56(n-1)*(n-2) + 21B*(n-1) +  6C,
        -27(n+0)*(n-1) - 6(B+c)*n  -   C,
          6(n+1)*(n+0) + 6c* (n+1)
    ]

    crit_flag = false
    for n ∈ 6:mord
        tmp = A(n)
        push!(α, -sum(tmp[1:7] .* α[end-6:end]) / last(tmp))

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

function transformed_z(a, b, c, z)
    zs = (;
        z = abs2(z),
        zoverzminus1a = abs2(z / (z - 1)),
        oneminusz = abs2(1 - z),
        oneminusoneoverz = abs2(1 - 1 / z),
        oneoverz = abs2(1 / z),
        oneoveroneminusz = abs2(1 / (1 - z))
    )

    tran = last(findmin(zs))

    if tran == :z
        if abs(c / a / b) < abs(c / (c - a) / (c - b))
            tran = :zalt
        end
    elseif tran == :zoverzminus1a
        if abs(c / a / (c - b)) < abs(c / (c - a) / b)
            tran = :zoverzminus1b
        end
    end

    return tran
end

function conformal_transformation_2f1(a, b, c, z; kwargs...)
    tran = transformed_z(a, b, c, z)

    return first(transformations[tran](a, b, c, z, conformal_2f1))
end
