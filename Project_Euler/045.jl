T(n) = n * ( n + 1) ÷ 2
P(n) = n * (3n - 1) ÷ 2
H(n) = n * (2n - 1)

function findtrianglenumber(n, N)
    for nT ∈ n:N
        t = T(nT)
        for nP ∈ 2:N
            p = P(nP)
            if p == t
                for nH ∈ 2:N
                    h = H(nH)
                    if h == t
                        return nT
                    end
                    if h > t
                        break
                    end
                end
            end
            if p > t
                break
            end
        end
    end

    return nothing
end

T(findtrianglenumber(286,10^5))
