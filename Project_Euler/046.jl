using Primes

function isgoldbach(n)
    if n % 2 == 0 || isprime(n)
        return true
    end

    sqrtn = floor(Int, sqrt(n / 2))

    for s ∈ 1 : sqrtn
        if isprime(n - 2s^2)
            return true
        end
    end

    return false
end

function goldbach(n)
    if n % 2 == 0 || isprime(n)
        return (n, 0)
    end

    sqrtn = floor(Int, sqrt(n / 2))

    for s ∈ 1 : sqrtn
        if isprime(n - 2s^2)
            return (n - 2s^2, s)
        end
    end

    return nothing
end

function get_nongoldbach(N)
    for n ∈ 9:2:N
        if !isgoldbach(n)
            return n
        end
    end

    return nothing
end
