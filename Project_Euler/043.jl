using Primes
using Combinatorics

function issubdivisible(d)
    for k = 2:8
        if parse(Int, join(d[k:k+2])) % prime(k - 1) != 0
            return false
        end
    end

    return true
end

function get_subdivisible_numbers()
    subdivisible_numbers = Vector{Int}()

    for n ∈ permutations(0:9)
        if issubdivisible(n)
            val = parse(Int, join(n))
            push!(subdivisible_numbers, val)
        end
    end

    return subdivisible_numbers
end

sum(get_subdivisible_numbers())
