sumDigits(N::Integer) =  sum(parse.(Int, collect(string(N))))

function powerSum(N::Integer)
    val = 2^big(N)

    return sumDigits(val)
end
