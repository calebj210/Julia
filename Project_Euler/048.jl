function bigpowersum(N)
    s = big(0)

    sum(big.(1:N).^(1:N))
end

bigpowersum(1000)
