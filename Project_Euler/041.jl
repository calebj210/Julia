using Primes

function ispandigital(n)
    d = sort(digits(n))

    return d == 1:length(d)
end

function getpandigitalprimes(N)
    ps = primes(N)

    pandigitalprimes = Vector{Int}()
    for p ∈ ps
        if ispandigital(p)
            push!(pandigitalprimes, p)
        end
    end

    return pandigitalprimes
end

getpandigitalprimes(987654321)
