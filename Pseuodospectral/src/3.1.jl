function fdweights(ξ, x, m)
    n = length(x) - 1
    # c = Array{Rational, 3}(undef, n + 1, n + 1, m + 1)
    c = zeros(Rational, n+1, n+1, m+1)
    c[1,1,1] = 1
    c1 = 1
    c4 = x[1] - ξ
    for i ∈ 1:n
        mn = min(i,m)
        c2 = 1
        c5 = c4
        c4 = x[i+1] - ξ
        for j ∈ 0:i-1
            c3 = x[i+1] - x[j+1]
            c2 *= c3
            if i < m
                c[j+1,i,i+1] = 0
            end
            c[j+1,i+1,1] = c4*c[j+1,i,1] // c3
            for k ∈ 1:mn
                c[j+1,i+1,k+1] = (c4*c[j+1,i,k+1] - k*c[j+1,i,k]) // c3
            end
        end
        c[i+1,i+1,1] = -c1*c5*c[i,i,1] // c2
        for k ∈ 1:mn
            c[i+1,i+1,k+1] = c1*(k*c[i,i,k] - c5*c[i,i,k+1]) // c2
        end
        c1 = c2
    end

    return c
end
