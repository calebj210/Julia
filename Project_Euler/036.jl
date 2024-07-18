function count_double_polindromes(N)
    palindromesum = 0
    for n ∈ 1 : N
        d10 = digits(n)
        d2  = digits(n, base = 2)
        if d10 == reverse(d10) && d2 == reverse(d2)
            palindromesum += n
        end
    end

    return palindromesum
end
