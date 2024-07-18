function ispandigital(n)
    d = sort(digits(n))

    return d == 1:9
end

function catprod(n, r)
    return parse(Int, join(n * (1:r)))
end

function get_largest_pandigital()
    pmax = 0
    amax = 0
    rmax = 0

    for r ∈ 8:-1:2
        for a ∈ 1:10^6
            n = catprod(a,r)
            if ispandigital(n) && n > pmax
                pmax = n
                amax = a
                rmax = r
            elseif n > 987654321
                break
            end
        end
    end

    return (amax, rmax, pmax)
end

get_largest_pandigital()
