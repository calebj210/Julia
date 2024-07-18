using Primes

function get_circular_primes(N::Int)
    circular_primes = Vector{Int}()

    for p ∈ primes(N)
        primeflag = true
        pdigits = reverse(digits(p))

        for n = 1:length(pdigits) - 1
            circp = parse(Int, join(circshift(pdigits, n)))
            if !isprime(circp)
                primeflag = false
                break
            end
        end

        if primeflag
            push!(circular_primes, p)
        end
    end

    return circular_primes
end

length(get_circular_primes(10^6))
