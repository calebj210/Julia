using Primes

function largestPrimeFactor(n::Integer)
    factors = factor(n)
    return maximum(keys(factors))
end
