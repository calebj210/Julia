using Primes

function evenlyDivisible(N::Integer)
    factors = [factor(n) for n ∈ 1 : N]

    result = reduce(factors; init = Dict()) do acc, dict
        mergewith(max, acc, dict)
    end

    return prodfactors(result)
end
