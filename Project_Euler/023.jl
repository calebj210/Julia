using Primes
using Base.Threads

function isabundant(n::Integer)
    divs = divisors(n)[1:end - 1]
    
    return sum(divs) > n
end

function isabundantsum(n::Integer)
    for k ∈ 12:n
        if !isabundant(k)
            continue                # Ignore non abundant numbers
        end

        if isabundant(n - k)
            return true
        end
    end
    
    return false
end

function nonabundant_sums()
    notasums = Vector{Integer}()
    for n = 1:28123
        if !isabundantsum(n)
            push!(notasums, n)
        end
    end

    return notasums
end

sum(nonabundant_sums())
