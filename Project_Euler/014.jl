using Base.Threads

function collatzLength(N::Integer)
    l = 1
    
    while N != 1
        N = iseven(N) ? N ÷ 2 : 3N + 1
        l += 1
    end

    return l
end

function maxLength(N::Integer)
    maxN = 1
    maxL = 1
    @threads for n ∈ 1 : N
        l = collatzLength(n)

        if l > maxL
            maxN = n
            maxL = l
        end
    end

    return (maxN, maxL)
end

maxLength(999999)
