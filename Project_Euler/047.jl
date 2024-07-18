using Primes

factorcount(n) = length(factor(n))

function first_conescutive_prime(n, N)
    for k ∈ 1 : N - n + 1
        if all(factorcount.(k:k + n - 1) .== n)
            return k
        end
    end

    return nothing
end
