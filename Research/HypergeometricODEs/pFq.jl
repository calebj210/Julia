#=
#   ODE approach for computing hypergeometric functions
# 
# Author: Caleb Jacobs
# DLM: January 30, 2025
=#

using Polynomials
using Polynomials.RationalFunctionFit: pade_fit

using MathLink, SparseArrays
using SpecialFunctions
using ArbNumerics: gamma, hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat
import SpecialFunctions.gamma

function sgn(x)
    iszero(x) ? -one(x) : sign(x)
end

gamma(z::Complex{BigFloat}) = gamma(ArbComplex(z))
poc(a, n) = iszero(n) ? 1 : prod(a + k for k ∈ 0:n - 1)

pade(p::AbstractPolynomial, m::Int, n::Int) = //(pade_fit(p, m, n)...)

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

johansson_2f1(a, b, c, z; bits = 106)::ComplexF64 = arb_2f1(ArbComplex.((a, b, c, z), bits = bits)...)

function maclaurin_2f1(a::Ta, b::Tb, c::Tc, z::Tz, N = 1000) where {Ta, Tb, Tc, Tz}
    coef_type = promote_type(Ta, Tb, Tc, Tz)

    S = coeff = one(coef_type)

    coeffs = Vector{coef_type}()
    push!(coeffs, coeff)

    zn = one(coef_type)

    # println("|c_0| = $(abs(coeff * zn))")
    for n ∈ 1 : N
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n
        push!(coeffs, coeff)

        zn *= z

        if abs(coeff*zn) <= eps(abs(S))
            break
        end

        S  += coeff * zn
        # println("|c_$(n)| = $(abs(coeff * zn))")
    end

    p  = Polynomial(coeffs)
    dp = derivative(p)

    return [p(z), dp(z)]
end

function recursive_2f1(a::Ta, b::Tb, c::Tc, z0::Tz, f0, h, N) where {Ta, Tb, Tc, Tz}
    Tcomplex = real(z0) >= 1 ? ComplexF64 : Float64             # Passed branch point 
    coeffs = Vector{promote_type(Ta, Tb, Tc, Tz, Tcomplex)}()
    push!(coeffs, f0...)

    # 10 flop optimization for 3-term recurrence
    a0 = -a * b; a1 =  1 - a - b; a2 = -2
    b0 = b1 = c - (1 + a + b) * z0; b2 = 2 - 4z0
    c0 = c1 = c2 = 2z0 * (z0 - 1)
    
    S = f0[1] + f0[2] * h
    hn = h
    for n = 3:N + 1
        # Compute next coefficient
        coeff = (a0 * coeffs[n - 2] + b0 * coeffs[n - 1]) / c0

        hn *= h
        criteria = abs(coeff * hn)
        if criteria <= eps(abs(S))
            push!(coeffs, coeff)
            break
        elseif isnan(criteria) || isinf(criteria)
            # @warn "Method diverged"
            break
        end

        push!(coeffs, coeff)
        S += coeff * hn

        # Update recurrence values
        a1 += a2; a0 += a1
        b1 += b2; b0 += b1
        c1 += c2; c0 += c1
    end 

    Tn  = Polynomial(coeffs)
    Tnp = derivative(Tn)

    return [Tn(h), Tnp(h)]
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

function taylor_terms(a, b, c, z0, p, h)
    coeffs = taylor_coefficients(a, b, c, z0, p)

    terms = coeffs .* h.^(0:(p - 1))

    fig = Figure()
    axt = Axis(fig[1,1], yscale = log10, ylabel = L"|a_n h^n|")
    axs = Axis(fig[2,1], yscale = log10, ylabel = L"|S_n|")
    axe = Axis(fig[3,1], yscale = log10, ylabel = L"\varepsilon(|S_n|)")
    scatter!(axt, abs.(terms))
    scatter!(axs, abs.(cumsum(terms)))
    scatter!(axe, eps.(abs.(cumsum(terms))))
    resize_to_layout!(fig)

    return (fig, terms)
end

function taylor_2f1(a, b, c, z::Number; H = Inf, N = 1000, order = 1000)
    if abs(z) <= .3
        return maclaurin_2f1(a, b, c, z, N)[1]
    end

    tmp = abs(c / (a * b))
    z0 = sign(z) * min(tmp, 0.3)
    # z0 = sign(z) * .3
    if real(z) > 1
        z0 += imag(z) > 0 ? 0.3im : -0.3im
    end
    dir = sign(z - z0)

    zn = z0
    fn = maclaurin_2f1(a, b, c, z0, N)

    n = 1
    while !isapprox(zn, z) && n < N
        r = abs(zn - 1)
        # h = dir * min(r / exp(2), abs(z - zn), H)           # Step size based on Jorba and Zou 2005
        h = dir * min(abs(1/fn[2]), r / exp(2), abs(z - zn), H)           # Step size based on Jorba and Zou 2005

        display(abs(fn[1]))
        display(h)
        fn = recursive_2f1(a, b, c, zn, fn, h, order + 1)   # Order increase so derivative hits the desired order
        display(zn)

        zn += h
        n += 1
    end

    return fn[1]
end

function _2f1(a, b, c, z::Number; H = 0.1, N = 1000, order = 1000)
    if real(z) <= 0.5 && abs(z) <= 1
        f = taylor_2f1(a, b, c, z, H = H, N = N, order = order)
    elseif abs(z) >= 1 && abs(z - 1) >= 1
        if isinteger(b - a)
            f = int_ab_2f1(a, b, c, z)
        else
            f1 = taylor_2f1(a, a - c + 1,  a - b + 1, 1 / z, H = H, N = N, order = order)
            f2 = taylor_2f1(b, b - c + 1, -a + b + 1, 1 / z, H = H, N = N, order = order)

            g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
            g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)

            f = g1 * f1 + g2 * f2
        end
    else
        if isinteger(c - a - b)
            f = int_abc_2f1(a, b, c, z)
        else
            f1 = taylor_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z, H = H, N = N, order = order)
            f2 = taylor_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z, H = H, N = N, order = order)

            g1 = gamma(c) / gamma(c - a) * gamma(c - a - b) / gamma(c - b) * z^(-a)
            g2 = gamma(c) / gamma(a)     * gamma(a + b - c) / gamma(b)     * (1 - z)^(c - a - b) * z^(a - c)

            f = g1 * f1 + g2 * f2
        end
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

function sparse_taylor_coefficients(a, b, c, z0, h, c0, c1, N)
    A(n) = -(n + a) * (n + b) * h^2
    B(n) =  (n + 1) * (n + c  - (1 + a + b + 2n) * z0) * h
    C(n) =  (n + 1) * (n + 2) * (1 - z0) * z0
    s0 = c0
    s1 = c1 * h

    As = A.(2:N - 2)
    Bs = B.(1:N - 2)
    Cs = C.(0:N - 2)
    
    LHS = spdiagm(-2 => As, -1 => Bs, 0 => Cs)
    RHS = sparsevec(Dict(1 => -A(0) * s0 - B(0) * s1, 2 => -A(1) * s1), N - 1)

    return [s0; s1; LHS\RHS]
end
