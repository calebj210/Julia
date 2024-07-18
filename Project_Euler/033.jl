function isdigitcancelling(num, den)
    n = digits(num)
    d = digits(den)

    if n[1] != d[2]
        return false
    end
    
    return (9n[2] + n[1]) * d[1] == 10n[2] * d[2] 
end

function find_digit_cancelling_fractions()
    digit_cancelling_fractions = Vector{Tuple{Int, Int}}()

    for den ∈ 11:99
        for num ∈ 10:den - 1
            if isdigitcancelling(num, den)
                push!(digit_cancelling_fractions, (num, den))
            end
        end
    end

    return digit_cancelling_fractions
end

function digitcancellingprod()
    fractions = find_digit_cancelling_fractions()
    fracs = [Rational(f...) for f ∈ fractions]

    return prod(fracs)
end

digitcancellingprod()
