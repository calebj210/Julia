using Base.Threads
using Primes

function recurring_digits_unit_fraction(denominator)
    # Error handling
    if denominator <= 0
        error("Denominator must be a positive integer.")
    end

    # Check if the fraction terminates (no recurring digits)
    if denominator == 1 || (denominator % 2 == 0 && denominator % 5 == 0)
        return 0
    end

    remainders = Int[]
    remainder = 1 % denominator  # Start with remainder 1

    while remainder ∉ remainders && remainder != 0
        push!(remainders, remainder)
        remainder = (remainder * 10) % denominator
    end

    return remainder == 0 ? 0 : length(remainders) - findfirst(isequal(remainder), remainders) + 1
end

function findMaxRecurringCount(N::Integer)
    maxCount = 0
    maxDenom = 0

    denoms = primes(N)

    @threads for d ∈ denoms
        count = recurring_digits_unit_fraction(d)
        if count > maxCount
            maxCount = count
            maxDenom = d
        end
    end

    return maxDenom
end
