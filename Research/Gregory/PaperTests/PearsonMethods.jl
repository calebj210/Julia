#=
# Implementation of hypergeometric methods discussed in 
#   Pearson, J.W., Olver, S. & Porter, M.A. Numerical methods for the computation of the confluent and Gauss hypergeometric functions. Numer Algor 74, 821–866 (2017).
# Author: Caleb Jacobs
# DLM: April 15, 2024
=#

using SpecialFunctions

"""
    2F1 Taylor series method (a)
"""
function taylorA(a, b, c, z, tol)
    if c == round(c, RoundToZero) && real(c) <= 0
        a1 = gamma(a - c + 1) / gamma(a) * gamma(b - c + 1) / gamma(b) / gamma(2 - c) * z^(1 - c)
        b1 = a1

        for j ∈ 1 : 500
            a1 = (a - c + j) * (b - c + j) / j * z / (j - c + 1) * a1
            b1 += a1
            
            # Stopping criteria
            if abs(a1) / abs(b1) < tol
                break
            end

            # No convergence
            if j == 500
                print("Taylor A did not converge in 500 terms")

                return
            end
        end
    else
        a1 = b1 = 1

        for j = 1 : 500
            a1 = (a + j - 1) * (b + j - 1) / (c + j - 1) * z / j * a1
            b1 += a1

            if abs(a1) / abs(b1) < tol
                break
            end

            if j == 500
                print("Taylor A did not converge in 500 terms")

                return
            end
        end
    end


    return b1
end

"""
    2F1 Taylor series method (b)
"""
function taylorB(a, b, c, z, tol)
    r = zeros(ComplexF64, 2)
    r[1] = a * b / c
    r[2] = (a + 1) * (b + 1) / 2 / (c + 1)

    A = zeros(ComplexF64, 2)
    A[1] = 1 + z * r[1]
    A[2] = A[1] + z^2 * a * b / c * r[2]

    for j = 3 : 500
        push!(r, (a + j - 1) * (b + j - 1) / j / (c + j - 1))

        push!(A, A[j - 1] + (A[j - 1] - A[j - 2])* r[j] * z)

        if abs(A[j] - A[j - 1]) / abs(A[j - 1]) < tol && abs(A[j - 1] - A[j - 2]) / abs(A[j - 2]) < tol
            break
        end

        if j == 500
            print("Taylor B did not converge in 500 terms")

            return
        end
    end

    return A[end]
end

"""
    2F1 Single fraction method
"""
function singleFraction(a, b, c, z, tol)
    a1 = ComplexF64.([0, c])
    b1 = ComplexF64.([1, a * b * z])
    c1 = ComplexF64.([1, c])
    d1 = ComplexF64.([1, (c + a * b * z) / c])
    
    for j = 3 : 500
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

        if j == 500
            print("Single fraction did not converge in 500 terms")
        end
    end

    return d1[end]
end
