#==
# Bisection method for root finding
# Author: Caleb Jacobs
# Date last modified: 24-08-2021
==#

f(x) = (x-1) * (x-3) * (x-5)

function bisect(f, a, b; maxIts = 10, tol = 10^(-5))
    mid = 0.0

    for i ∈ 1:maxIts
        mid = (a + b) / 2
        
        if f(mid) == 0 || (b - a) / 2 < tol
            return mid
        end

        if sign(f(mid)) == sign(f(a))
            a = mid
        else
            b = mid
        end
    end

    return mid
end

display(bisect(f, 0.0, 20.0, maxIts = 100, tol = 10^(-10)))
