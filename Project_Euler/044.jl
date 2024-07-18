P(n) = n * (3n - 1) ÷ 2

function ispentagonal(n)
    return isinteger((sqrt(24n + 1) + 1) / 6)
end

function pentagonaldifference(N)
    for k ∈ 1 : N
        for n ∈ 1 : k
            if ispentagonal(P(n) + P(k)) && ispentagonal(abs(P(n) - P(k)))
                return (k, n, abs(P(n) - P(k)))
            end
        end
    end
    
    return nothing
end

pentagonaldifference(10^4)
