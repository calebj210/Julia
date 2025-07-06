#=
# Implementation of hypergeometric methods discussed in 
#   Pearson, J.W., Olver, S. & Porter, M.A. Numerical methods for the computation 
#   of the confluent and Gauss hypergeometric functions. Numer Algor 74, 821–866 (2017).
#
# Author: Caleb Jacobs
# DLM: October 2, 2024
=#

using SpecialFunctions
using FastGaussQuadrature

"2F1 Taylor series method (a)"
function taylor_a_2f1(a, b, c, z, tol = 1e-15)
    if c == round(c, RoundToZero) && real(c) <= 0
        a1 = gamma(a - c + 1) / gamma(a) * gamma(b - c + 1) / gamma(b) / gamma(2 - c) * z^(1 - c)
        b1 = a1

        for j ∈ 1 : 1000
            a1 = (a - c + j) * (b - c + j) / j * z / (j - c + 1) * a1
            b1 += a1
            
            # Stopping criteria
            if abs(a1) / abs(b1) < tol
                break
            end

            # No convergence
            if j == 1000
#                 print("Taylor A did not converge in 1000 terms")

                return complex(NaN)
            end
        end
    else
        a1 = b1 = 1

        for j = 1 : 1000
            a1 = (a + j - 1) * (b + j - 1) / (c + j - 1) * z / j * a1
            b1 += a1

            if abs(a1) / abs(b1) < tol
                break
            end

            if j == 1000
#                 print("Taylor A did not converge in 1000 terms")

                return complex(NaN)
            end
        end
    end


    return b1
end

"2F1 Taylor series method (b)"
function taylor_b_2f1(a, b, c, z, tol = 1e-15)
    r = zeros(ComplexF64, 2)
    r[1] = a * b / c
    r[2] = (a + 1) * (b + 1) / 2 / (c + 1)

    A = zeros(ComplexF64, 2)
    A[1] = 1 + z * r[1]
    A[2] = A[1] + z^2 * a * b / c * r[2]

    for j = 3 : 1000
        push!(r, (a + j - 1) * (b + j - 1) / j / (c + j - 1))

        push!(A, A[j - 1] + (A[j - 1] - A[j - 2])* r[j] * z)

        if abs(A[j] - A[j - 1]) / abs(A[j - 1]) < tol && abs(A[j - 1] - A[j - 2]) / abs(A[j - 2]) < tol
            break
        end

        if j == 1000
#             print("Taylor B did not converge in 1000 terms")

            return complex(NaN)
        end
    end

    return A[end]
end

" 2F1 Single fraction method "
function single_fraction_2f1(a, b, c, z, tol = 1e-15)
    a1 = ComplexF64.([0, c])
    b1 = ComplexF64.([1, a * b * z])
    c1 = ComplexF64.([1, c])
    d1 = ComplexF64.([1, (c + a * b * z) / c])
    
    for j = 3 : 1000
        push!(a1, (a1[j - 1] + b1[j - 1]) * (j - 1) * (c + j - 2))
        push!(b1, b1[j - 1] * (a + j - 2) * (b + j - 2) * z)
        push!(c1, c1[j - 1] * (j - 1) * (c + j - 2))

        # Stop if any of these are infinitely large
        if isinf(a1[j]) || isinf(b1[j]) || isinf(c1[j])
            break
        end

        push!(d1, (a1[j] + b1[j]) / c1[j])

        if abs(d1[j] - d1[j - 1]) / abs(d1[j - 1]) < tol && abs(d1[j - 1] - d1[j - 2]) / abs(d1[j - 2]) < tol
            break
        end

        if j == 1000
#             print("Single fraction did not converge in 1000 terms")

            return complex(NaN)
        end
    end

    return d1[end]
end

" 2F1 Buhring's expansion formula "
function buhring_2f1(a, b, c, z, z0 = 0.0 + 0im, tol = 1e-15)
    d = zeros(ComplexF64, 2)
    d[1] = (1 + a - b)^(-1) * a * ((a + 1) * (1 - 2z0) + (a + b + 1) * z0 - c)
    d[2] = (2 * (2 + a - b))^(-1) * (a + 1) * (((a + 2) * (1 - 2z0) + (a + b + 1) * z0 - c) *
           (1 + a - b)^(-1) * a * ((a + 1) * (1 - 2z0) + (a + b + 1) * z0 - c) + z0 * (1 - z0) * a)

    a1 = 1 + d[1] * (z - z0)^(-1) + d[2] * (z - z0)^(-2)
    b1 = a1
    for n ∈ 3 : 1000
        push!(d, (n * (n + a - b))^(-1) * (n + a - 1) *
               (((n + a) * (1 - 2z0) + (a + b + 1) * z0 - c) * d[n - 1] + z0 * (1 - z0) * (n + a - 2) * d[n - 2]))
        a1 = d[n] * (z - z0)^(-n)
        b1 += a1

        if abs(a1) / abs(b1) < tol
            break
        end

        if n == 1000
#             print("Buhring did not converge to tolerance in 1000 iterations")
            
            return complex(NaN)
        end
    end

    e = zeros(ComplexF64, 2)
    e[1] = (1 - a + b)^(-1) * b * ((b + 1) * (1 - 2z0) + (a + b + 1) * z0 - c)
    e[2] = (2 * (2 - a + b))^(-1) * (b + 1) * (((b + 2) * (1 - 2z0) + (a + b + 1) * z0 - c) *
           (1 - a + b)^(-1) * b * ((b + 1) * (1 - 2z0) + (a + b + 1) * z0 - c) + z0 * (1 - z0) * b)

    a2 = 1 + e[1] * (z - z0)^(-1) + e[2] * (z - z0)^(-2)
    b2 = a2
    for n ∈ 3 : 1000
        push!(e, (n * (n - a + b))^(-1) * (n + b - 1) * 
                 (((n + b) * (1 - 2z0) + (a + b + 1) * z0 - c) * e[n - 1] + z0 * (1 - z0) * (n + b - 2) * e[n - 2]))
        a2 = e[n] * (z - z0)^(-n)
        b2 += a2

        if abs(a2) / abs(b2) < tol
            break
        end

        if n == 1000
#             print("Buhring did not converge to tolerance in 1000 iterations")

            return complex(NaN)
        end
    end

    if real(b) < 0 && isinteger(b) || real(a) < 0 && isinteger(a) || real(c) < 0 && isinteger(c) || isinteger(a - b)
        return NaN
    else
        return gamma(c) * (gamma(b - a) / gamma(b) / gamma(c - a) * (z0 - z)^(-a) * b1 +
                           gamma(a - b) / gamma(a) / gamma(c - b) * (z0 - z)^(-b) * b2)
    end
end

" 2F1 Gauss-Jacobi quadrature "
function gauss_jacobi_quadrature_2f1(a, b, c, z, N = 150)
    x = Vector{ComplexF64}(undef, N)
    w = Vector{ComplexF64}(undef, N)

    try
        x, w = gaussjacobi(N, c - b - 1, b - 1)
    catch
        return NaN
    end

    return [gamma(c) / gamma(b) / gamma(c - b) / (2^(c - 1)) *
           sum(w .* (1 .- z / 2 * (x .+ 1)).^(-a)) for z ∈ z]
end
gauss_jacobi_quadrature_2f1(a, b, c, z::Number, N = 150) = gauss_jacobi_quadrature_2f1(a, b, c, [z])[1]
