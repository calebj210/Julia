using Primes

function isleftrightprime(n)
    d = reverse(digits(n))

    primeflag = true

    for k = 1:length(d)
        lr = parse(Int, join(d[1:k]))
        rl = parse(Int, join(d[end-k+1:end]))

        if !isprime(lr) || !isprime(rl)
            primeflag = false
            break
        end
    end

    return primeflag
end

function getleftrightprimes(N)
    leftrightprimes = Vector{Int}()

    for n = primes(10, N)
        if isleftrightprime(n)
            push!(leftrightprimes, n)
        end

        if length(leftrightprimes) == 11
            break
        end
    end

    return leftrightprimes
end

sum(getleftrightprimes(10^6))
