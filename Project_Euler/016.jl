function powerdigitsum(n)
    p = big(2)^n

    return sum(digits(p))
end

powerdigitsum(1000)
