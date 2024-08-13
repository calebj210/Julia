using Primes
using Combinatorics

function ispermutation(a, b)
    da = sort(digits(a))
    db = sort(digits(b))

    return da == db
end

function isprimepermutation(n)
    for pd ∈ permutations(reverse(digits(n)))
        p = parse(Int, join(pd))

        if p <= n || !isprime(p)
            continue
        end

        q = 2p - n
        if ispermutation(n, q) && isprime(q)
            return true
        end
    end

    return false
end

function getprimepermutation(n)
    for pd ∈ permutations(reverse(digits(n)))
        p = parse(Int, join(pd))

        if p <= n || !isprime(p)
            continue
        end

        q = 2p - n
        if ispermutation(n, q) && isprime(q)
            return [n,p,q]
        end
    end

    return false
end

function getprimeperms()
    ps = primes(10^3, 10^4)

    primeperms = Vector{Int}()

    for p ∈ ps
        if length(p) >= 2
            break
        end
        if isprimepermutation(p)
            push!(primeperms, p)
        end
    end

    return primeperms
end

function readprimeperms()
    ps = getprimeperms()

    for p ∈ ps
        println(join(getprimepermutation(p)))
    end
end

readprimeperms()
