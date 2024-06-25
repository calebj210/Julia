using Primes

function primeSum(N::Integer)
    ps = primes(2, N - 1)

    return sum(ps)
end
