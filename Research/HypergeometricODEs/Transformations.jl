#=
#   Transformations for 2F1
#
# Author: Caleb Jacobs
# DLM: June 5, 2025
=#

using SpecialFunctions

"""
    int_abc_2f1(a, b, c, z; N = 150)
Evaluate 2F1 when `c`-`a`-`b` is an integer using the integer limiting case formula from https://dlmf.nist.gov/15.8.E11
"""
function int_abc_2f1(a, b, c, z; N = 100, tol = 1e-15)
    m = c - a - b                           # Get integer parameter
    if !isinteger(m)
        error("c - a - b should be an integer")
    end

    negflag = false
    if m < 0
        m *= -1
        a -= m
        b -= m
        
        negflag = true
    end

    coef = gamma(m) / gamma(b + m)
    sum1 = coef
    for k ∈ 1:m - 1
        coef *= (a + k - 1) / (m - k) / k * (b + m - k) * (1 - 1 / z)
        sum1 += coef

        if abs(coef / sum1) <= tol
            break
        end
    end
    sum1 *= z^(-a) / gamma(a + m)

    coef = 1 / gamma(m + 1) / gamma(b) * (1 - 1 / z)^m
    prd  = log((1 - z) / z) - digamma(1) - digamma(m + 1) + digamma(a + m) + digamma(b)
    sum2 = coef * prd
    for k ∈ 1:N
        coef *= (a + m + k - 1) / k / (k + m) * (b - k) * (-1) * (1 - 1 / z)
        prd = log((1 - z) / z) - digamma(k + 1) - digamma(k + m + 1) + digamma(a + k + m) + digamma(b - k) 

        sum2 += coef * prd

        if abs(coef * prd / sum2) <= tol
            break
        end
    end
    sum2 *= -z^(-a) / gamma(a)

    if negflag
        return gamma(c) * (1 - z)^(-m) * (sum1 + sum2)
    else
        return gamma(c) * (sum1 + sum2)
    end
end

"""
    int_ab_2f1(a, b, c, z; N = 150)
Evaluate 2F1 when `b`-`a` is an integer using the integer limiting case formula from https://dlmf.nist.gov/15.8.E8
"""
function int_ab_2f1(a, b, c, z; N = 100, tol = 1e-15)
    m = b - a                               # Get integer parameter
    if !isinteger(m)
        error("b - a should be an integer")
    end

    if m < 0
        m *= -1
        tmp = a
        a = b
        b = tmp
    end

    coef = gamma(m) / gamma(c - a)
    sum1 = coef
    for k ∈ 1:m - 1
        coef *= (a + k - 1) / (m - k) / k * (c - a - k) * (1 / z)
        sum1 += coef

        if abs(coef / sum1) <= tol
            break
        end
    end
    sum1 *= (-z)^(-a) / gamma(a + m)

    coef = 1 / gamma(m + 1) / gamma(c - a - m) * (1 / z)^m
    prd  = log(-z) + digamma(1) + digamma(m + 1) - digamma(a + m) - digamma(c - a - m)
    sum2 = coef * prd
    for k ∈ 1:N
        coef *= (a + m + k - 1) / k / (k + m) * (c - a - k - m) * (-1) * (1 / z)
        prd = log(-z) + digamma(k + 1) + digamma(k + m + 1) - digamma(a + k + m) - digamma(c - a - k - m) 

        sum2 += coef * prd

        if abs(coef * prd / sum2) <= tol
            break
        end
    end
    sum2 *= (-z)^(-a) / gamma(a)

    return gamma(c) * (sum1 + sum2)
end

function maclaurin_2f1(a, b, c, z, N = 1000; tol = eps() / 2)
    S = zn = coeff = 1.0
    dS = 0.0

    crit_flag = false
    for n ∈ 1:N
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        dS += coeff * zn * n
        zn *= z
        val = coeff * zn
        if abs(val / S) <= tol / 2
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

    return [S, dS]
end

function zalt_2f1(a, b, c, z)
    f, df = maclaurin_2f1(c - a, c - b, c, z)

    g = (1 - z)^(c - a - b)

    return [
        g * f,
        g * ((a + b - c) * f / (1 - z) + df)
    ]
end

function zoverzminus1a_2f1(a, b, c, z)
    f, df = maclaurin_2f1(a, c - b, c, z / (z - 1))

    g = (1 - z)^(-a)

    return [
        g * f,
        g * (a * f - df / (1 - z)) / (1 - z)
    ]
end

function zoverzminus1b_2f1(a, b, c, z)
    f, df = maclaurin_2f1(c - a, b, c, z / (z - 1))

    g = (1 - z)^(-b)

    return [
        g * f,
        g * (b * f - df / (1 - z)) / (1 - z)
    ]
end

function oneminusoneoverz_2f1(a, b, c, z)
    f1, df1 = maclaurin_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z)
    f2, df2 = maclaurin_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z)

    g1 = gamma(c - a - b) * gamma(c) / gamma(c - a) / gamma(c - b) * z^(-a)
    g2 = gamma(a + b - c) * gamma(c) / gamma(a) / gamma(b) * (1 - z)^(c - a - b) * z^(a - c)

    return [
        g1 * f1 + g2 * f2,
        g1 * (-a * f1 + df1 / z) / z +
        g2 * ((a + b - c) * f2 / (1 - z) + ((a - c) * f2 + df2 / z) / z)
    ]
end

function oneminusz_2f1(a, b, c, z)
    f1, df1 = maclaurin_2f1(a, b,         a + b - c + 1, 1 - z)
    f2, df2 = maclaurin_2f1(c - a, c - b, c - a - b + 1, 1 - z)
    
    g1 = gamma(c - a - b) * gamma(c) / gamma(c - a) / gamma(c - b)
    g2 = gamma(a + b - c) * gamma(c) / gamma(a) / gamma(b) * (1 - z)^(c - a - b)

    return [
        g1 * f1 + g2 * f2,
        -g1 * df1 + 
        g2 * ((a + b - c) / (1 - z) * f2 - df2)
    ]
end

function oneoverz_2f1(a, b, c, z)
    f1, df1 = maclaurin_2f1(a, a - c + 1, a - b + 1, 1 / z)
    f2, df2 = maclaurin_2f1(b, b - c + 1, b - a + 1, 1 / z)

    g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
    g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)

    return [
        g1 * f1 + g2 * f2, 
        (g1 * (-a * f1 - df1 / z) + 
         g2 * (-b * f2 - df2 / z)) / z
    ]
end

function oneoveroneminusz_2f1(a, b, c, z)
    f1, df1 = maclaurin_2f1(a, c - b, a - b + 1, 1 / (1 - z))
    f2, df2 = maclaurin_2f1(b, c - a, b - a + 1, 1 / (1 - z))

    g1 = gamma(b - a) * gamma(c) / gamma(b) / gamma(c - a) * (1 - z)^(-a)
    g2 = gamma(a - b) * gamma(c) / gamma(a) / gamma(c - b) * (1 - z)^(-b)

    return [
        g1 * f1 + g2 * f2,
        (g1 * (a * f1 + df1 / (1 - z)) +
         g2 * (b * f2 + df2 / (1 - z))) / (1 - z)
    ]
end

const transformations = (;
    z = maclaurin_2f1,
    zalt = zalt_2f1,

    zoverzminus1a = zoverzminus1a_2f1,
    zoverzminus1b = zoverzminus1b_2f1,

    oneminusz = oneminusz_2f1,
    oneminusoneoverz = oneminusoneoverz_2f1,

    oneoverz = oneoverz_2f1,
    oneoveroneminusz = oneoveroneminusz_2f1,
)
