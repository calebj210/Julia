function quadraticPrimeCount(a::Integer, b::Integer)
    n = 0
    while isprime(n^2 + a * n + b)
        n += 1
    end

    return n
end

function findBestQuadraticCoefficients(N::Integer)
    bestA = 0
    bestB = 0
    bestCount = 0

    for b ∈ primes(0, N)
        for a ∈ -(N - 1) : N - 1
            count = quadraticPrimeCount(a, b)
            if count > bestCount
                bestCount = count
                bestA = a
                bestB = b
            end
        end
    end

    return (bestA, bestB)
end
