using Primes

function findmaxprimesum(N)
    ps = primes(N)

    maxl = 1
    maxp = 2

    for n ∈ 1:length(ps)
        for k ∈ 1:length(ps) - n + 1
            s = sum(ps[k:k + n - 1])
            if s > N
                break
            end

            if isprime(s) && n > maxl
                maxp = s
                maxl = n
                break
            end
        end
    end

    return (maxp, maxl)
end
