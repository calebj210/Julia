function pythagorean_triples_by_perimeter(max_perimeter)
    triples = Vector{Tuple{Int,Int,Int}}()
    m = 2

    while m * (m + 1) <= max_perimeter  # New condition based on perimeter
        for n in 1:(m - 1)
            if (m + n) % 2 == 1 && gcd(m, n) == 1
                a, b, c = m^2 - n^2, 2 * m * n, m^2 + n^2

                if a + b + c <= max_perimeter  # Check perimeter
                    push!(triples, (a, b, c))

                    k = 2
                    while k * (a + b + c) <= max_perimeter
                        push!(triples, (k * a, k * b, k * c))
                        k += 1
                    end
                end
            end
        end
        m += 1
    end

    return triples
end

function countperimeters(triples, p)
    return count(sum.(triples) .== p)
end

function maxperimetersolutions(N)
    triples = pythagorean_triples_by_perimeter(N)

    mp = 0
    mc = 0
    for n ∈ 1 : N
        c = countperimeters(triples, n)
        if c > mc
            mp = n
            mc = c
        end
    end

    return (mp, mc)
end
