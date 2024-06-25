divisorSum(n::Integer) = sum(divisors(n)[1 : end - 1])

function amicableSum(N::Integer)
    divs = Dict{Int, Vector{Int}}()

    s = 0
    for n ∈ 2 : N
        d = divisorSum(n)
        if d >= N || d == n
            continue
        end

        if divisorSum(d) == n
            s += n
        end
    end

    return s
end

amicableSum(10000)
