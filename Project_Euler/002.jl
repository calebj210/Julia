function evenFibonacci(N::Integer)
    fib = [1, 2]
    s = 0

    while fib[2] <= N
        if fib[2] % 2 == 0
            s += fib[2]
        end

        tmp    = fib[2]
        fib[2] = sum(fib)
        fib[1] = tmp
    end

    return s
end
