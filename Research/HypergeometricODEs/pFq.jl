#=
#   ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: March 4, 2025
=#

using MathLink
using SpecialFunctions
using ArbNumerics: gamma, hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat
import SpecialFunctions.gamma

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

function maclaurin_2f1(a, b, c, z, N = 1000; tol = eps())
    S = zn = coeff = 1.0
    dS = 0.0

    for n ∈ 1:N
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        dS += coeff * zn * n
        zn *= z
        val = coeff * zn
        if abs(val / S) <= tol
        # if abs(val) <= eps(abs(S))
            break
        end
        S += val
    end

    return [S, dS]
end

function recursive_2f1(a, b, c, z0, f0, h, N; tol = eps())
    (c0, c1) = f0
    c1 *= h

    # 10 flop optimization for 3-term recurrence
    A0 = -a * b * h^2; A1 = (1 - a - b) * h^2; A2 = -2h^2
    B0 = B1 = (c - (1 + a + b) * z0) * h; B2 = (2 - 4z0) * h
    C0 = C1 = C2 = 2z0 * (z0 - 1)
    
    S  = c0 + c1
    dS = f0[2]
    for n = 3:N + 1
        # Compute next coefficient
        coeff = (A0 * c0 + B0 * c1) / C0

        dS += coeff / h * (n - 1)
        # criteria = abs(val)
        # if criteria <= eps(abs(S))
        criteria = abs(coeff / S)
        if criteria <= tol
            break
        elseif isnan(criteria) || isinf(criteria)
            # @warn "Method diverged"
            break
        end
        S += coeff

        c0 = c1
        c1 = coeff

        # Update recurrence values
        A1 += A2; A0 += A1
        B1 += B2; B0 += B1
        C1 += C2; C0 += C1
    end 

    return [S, dS]
end

function taylor_2f1(a, b, c, z::Number; N = 1000, order = 1000, step_max = Inf, init_max = exp(-1))
    # if abs(a * b) > 5 || abs(c) < 1
        init_max = min(abs(c / (a * b)), init_max)
    # end
    if abs(z) <= init_max
        return maclaurin_2f1(a, b, c, z, N)[1]
    end

    if real(z) > 1
        # z0 = init_max * im * (imag(z) > 0 ? 1 : -1)
        # znew = im * (imag(z) > 0 ? 1 : -1)
        znew = 1 + 3(imag(z) > 0 ? im : -im)
        z0 = init_max * sign(znew)
        # z0 = init_max * im * (imag(z) > 0 ? 1 : -1)
    else
        z0 = sign(z) * init_max
    end

    fn = maclaurin_2f1(a, b, c, z0, N)
    if real(z) > 1
        fn = taylor_init(a, b, c, z0, znew, fn)
        z0 = znew
    end
    dir = sign(z - z0)

    n = 1
    while !isapprox(z0, z) && n < N
        h_opt = abs(z0 - 1) / exp(2)
        # h_opt = abs(z0 - 1)
        h_end = abs(z0 - z)
        h_ord = abs(z0) / exp(2)
        # h_ord = sqrt(abs(2z0 * (1 - z0) / (a * b)))                    # sqrt(C / A)
        # h_ord = Inf
        # h_ord = abs(2z0 * (1 - z0) / (c - (1 + a + b) * z0))                  # C
        # h_ord = abs(2z0 * (1 - z0) * (c - (1 + a + b) * z0) / (a * b))        # C B / A

        h = dir * min(h_opt, h_end, h_ord, step_max)        # Step size based on Jorba and Zou 2005

        fn = recursive_2f1(a, b, c, z0, fn, h, order + 1)   # Order increase so derivative hits the desired order

        z0 += h
        n += 1
    end

    return fn[1]
end

function taylor_init(a, b, c, z0, z, f; max_step_size = Inf, max_steps = 1000, max_order = 1000)
    dir = sign(z - z0)
    fn = f
    n = 1
    while !isapprox(z0, z) && n < max_steps
        h_opt = abs(z0 - 1) / exp(2)
        h_end = abs(z0 - z)
        h_ord = abs(z0) / exp(2)
        # h_ord = sqrt(abs(2z0 * (1 - z0) / (a * b)))             # sqrt(C / A)
        # h_ord = Inf

        h = dir * min(h_opt, h_end, h_ord, max_step_size)       # Step size based on Jorba and Zou 2005

        fn = recursive_2f1(a, b, c, z0, fn, h, max_order + 1)   # Order increase so derivative hits the desired order

        z0 += h
        n += 1
    end

    return fn
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
    elseif  abs(z) >= 3 && abs(1 - z) >= 3
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

function taylor_coefficients(a,b,c,z0,p)
    F, dF = maclaurin_2f1(a,b,c,z0)

    A0 = a + b - a * b - 1; A1 =  3 - a - b; A2 = -2    # Initial A(n)
    B0 = 0; B1 = c - (a + b - 3) * z0 - 2; B2 = 2 - 4z0 # Initial B(n)
    C0 = C1 = 0; C2 = 2z0 * (z0 - 1)                    # Initial C(n)
    an = [F, dF]                                        # a_0 and a_1
    
    # 10-flop recursion
    for n = 3:p
        A1 += A2; A0 += A1                              # Update A(n)
        B1 += B2; B0 += B1                              # Update B(n)
        C1 += C2; C0 += C1                              # Update C(n)
        push!(an, (A0 * an[n - 2]  +                    # Next a_n
                   B0 * an[n - 1]) / C0)
    end
    return an
end

function As_Bs_and_Cs(a, b, c, z0, N = 50)
    A(n) = -(a+n)*(b+n)
    B(n) = (1+n)*(c-(1+a+b+2n)*z0+n)
    C(n) = (1+n)*(2+n)*(1-z0)*z0

    abcs = [A.(0:N) B.(0:N) C.(0:N)]

    return abcs
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
