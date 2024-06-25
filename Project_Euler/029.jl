function distinct_powers(N::Integer)
    a = big.(2:N)

    return length(union(vec(big.(a) .^ b')))
end

distinct_powers(100)
