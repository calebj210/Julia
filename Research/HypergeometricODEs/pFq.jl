#=
# ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: April 15, 2025
=#

using MathLink
using SpecialFunctions
using LinearAlgebra
using ArbNumerics: gamma, hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat
import SpecialFunctions.gamma
using Polynomials

function sgn(x)
    iszero(x) ? -one(x) : sign(x)
end

mathematica_2f1(a, b, c, z) = 
    try 
        Complex(weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z).args...)
    catch 
        Real(weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z))
    end
    
mathematica_pfq(a, b, z) = 
    try 
        Complex(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z).args...)
    catch 
        Real(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z))
    end

johansson_2f1(a, b, c, z; bits = 512)::ComplexF64 = arb_2f1(ArbComplex.((a, b, c, z), bits = bits)...)

function get_direction(z0, z, f, df, branch = false)
    if branch && real(z) > 1 && real(z0) < 1 
        h_str = angle(1 + 2sgn(imag(z)) * im - z0)  # Go to 1 ± i first when navigating the branch point
    else
        h_str = angle(z - z0)                       # Straight path direction
        branch = false
    end
    h_arg = angle(im * f / df)                           # Constant phase direction
    if h_arg < 0
        h_arg += π
    end

    h_dif = argmin(h -> abs2(h), h_arg - h_str .+ (-2:1) * π)

    return (cis(h_str + .1h_dif), branch)
end

function correct_taylor(a, b, c, f0, f1, h, z1; order = 1000, tol = eps(), max_iter = 5)
    f = recursive_2f1(a, b, c, z1, f1, -h, order)
    f10 = recursive_2f1(a, b, c, z1, [1.0, 0.0im], -h, order)
    f01 = recursive_2f1(a, b, c, z1, [0.0im, 1.0], -h, order)
    
    J = [f10(-h) f01(-h); derivative(f10)(-h) derivative(f01)(-h)]
    for _ ∈ 1:max_iter
        F = [f(-h); derivative(f)(-h)] - f0
        f1 -= J \ F
        if norm(F) / norm(f0) <= tol / 2
            break
        end
        f = recursive_2f1(a, b, c, z1, f1, -h, order)
    end

    return f1
end

function maclaurin_2f1(a, b, c, z, N = 1000; tol = eps())
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

function recursive_2f1(a, b, c, z0, f0, h, N; tol = eps())
    p = Polynomial(f0)
    sizehint!(p.coeffs, 50, shrink = false)

    (c0, c1) = f0

    # 10 flop optimization for 3-term recurrence
    A0 = -a * b; A1 = (1 - a - b); A2 = -2
    B0 = B1 = (c - (1 + a + b) * z0); B2 = (2 - 4z0)
    C0 = C1 = C2 = 2z0 * (z0 - 1)
    
    hn = h
    S = c0 + c1 * h
    crit_flag = false
    for n = 3:N + 1
        # Compute next coefficient
        coeff = (A0 * c0 + B0 * c1) / C0

        # criteria = abs(val)
        # if criteria <= eps(abs(S))
        hn *= h
        criteria = abs(coeff * hn / S)
        if criteria <= tol / 2
            if crit_flag
                break
            else
                crit_flag = true
            end
        elseif isnan(criteria) || isinf(criteria)
            # @warn "Method diverged"
            break
        elseif crit_flag
            crit_flag = false
        end
        S += coeff * hn

        c0 = c1
        c1 = coeff
        push!(p.coeffs, coeff)

        # Update recurrence values
        A1 += A2; A0 += A1
        B1 += B2; B0 += B1
        C1 += C2; C0 += C1
    end 

    return p
end

function taylor_2f1(a, b, c, z::Number; N = 1000, order = 1000, step_max = Inf, init_max = .3, backward = false)
    init_max = min(3abs(c / (a * b)), init_max)
    if abs(z) <= init_max
        return maclaurin_2f1(a, b, c, z, N)[1]
    end

    z0 = init_max * sign(z)

    fn = maclaurin_2f1(a, b, c, z0, N)

    branch = true
    for _ ∈ 1:N
        h_opt = abs(z0 - 1) * exp(-2)
        h_end = abs(z0 - z)

        h_ord = abs(z0) * exp(-2)
        # h_ord = Inf

        # h_rat = abs(fn[1] / fn[2])
        h_rat = Inf
        
        dir, branch = get_direction(z0, z, fn..., branch)
        step_size = min(h_opt, h_end, h_ord, h_rat, step_max)

        h = dir * step_size
        p = recursive_2f1(a, b, c, z0, fn, h, order)
        dp = derivative(p)
        if step_size !== h_end
            if !backward
                fn = [p(h), dp(h)]
                z0 += h
            else
                z0 += h
                fn = correct_taylor(a, b, c, fn, [p(h), dp(h)], h, z0, max_iter = 20)
            end
        else
            h = z - z0
            if backward
                fn = correct_taylor(a, b, c, fn, [p(h), dp(h)], h, z, max_iter = 20)
                return fn[1]
            else
                return p(h)
            end
        end
        
    end

    return fn[1]
end

function _2f1(a, b, c, z::Number; step_max = Inf, N = 1000, order = 1000)
    if real(z) >= 0.5 && abs(1 - z) <= 1
        if isinteger(c - a - b)
            f = int_abc_2f1(a, b, c, z)
            f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
        else
            f1 = taylor_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z, step_max = step_max, N = N, order = order)
            f2 = taylor_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z, step_max = step_max, N = N, order = order)

            g1 = gamma(c) / gamma(c - a) * gamma(c - a - b) / gamma(c - b) * z^(-a)
            g2 = gamma(c) / gamma(a)     * gamma(a + b - c) / gamma(b)     * (1 - z)^(c - a - b) * z^(a - c)

            f = g1 * f1 + g2 * f2
        end
    elseif abs(z) >= 1 && abs(1 - z) >= 1
    # if abs(z) >= 1 && abs(1 - z) >= 1
        if isinteger(b - a)
            f = int_ab_2f1(a, b, c, z)
            f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
        else
            f1 = taylor_2f1(a, a - c + 1,  a - b + 1, 1 / z, step_max = step_max, N = N, order = order)
            f2 = taylor_2f1(b, b - c + 1, -a + b + 1, 1 / z, step_max = step_max, N = N, order = order)

            g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
            g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)

            f = g1 * f1 + g2 * f2
        end
    else
        f = taylor_2f1(a, b, c, z, step_max = step_max, N = N, order = order)
    end

    return f
end

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
