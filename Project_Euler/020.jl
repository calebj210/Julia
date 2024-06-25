function factorial_sum(n::Integer)
    return sum(parse.(Integer, collect(string(factorial(big(n))))))
end

factorial_sum(100)
