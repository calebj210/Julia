function isdigitpowersum(n, p)
    powersum = sum(digits(n).^p)

    return n == powersum
end

function powersums(p, N)
    sums = Vector{Integer}()
    for n ∈ 2 : N 
        if isdigitpowersum(n, p)
            push!(sums, n)
        end
    end

    return sums
end

sum(powersums(5, 10^6))
