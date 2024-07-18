function getproduct(N)
    d = join(1:10^N)

    vals = parse.(Int, collect(d[10 .^(0:N)]))
    
    return prod(vals)
end
