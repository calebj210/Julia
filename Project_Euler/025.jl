function fibCount(N::Integer)
    f = big.([1, 1])
    n = 2
    
    while f[2] < big(10)^N
        tmp   = f[2]
        f[2] += f[1]
        f[1]  = tmp

        n += 1
    end

    return n
end
