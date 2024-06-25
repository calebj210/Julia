using Primes

function divisibleTriangle(N::Integer)
    tri(n) = n * (n + 1) ÷ 2

    factorLength = 0
    i = 1
    while factorLength <= N
        i += 1
        factorLength = prod(values(factor(tri(i))) .+ 1)
    end

    return tri(i)
end
