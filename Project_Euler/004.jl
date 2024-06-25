function largest_palindrome_product()
    largest_palindrome = 0
    for i in 999:-1:100  # Iterate from 999 down to 100
        for j in 999:-1:i  # Start j at i to avoid duplicate checks
            product = i * j
            if is_palindrome(product) && product > largest_palindrome
                largest_palindrome = product
            end
        end
    end
    return largest_palindrome
end

function is_palindrome(n)
    n_str = string(n)
    return n_str == reverse(n_str)
end

largest_palindrome = largest_palindrome_product()
println("The largest palindrome made from two 3-digit numbers is: ", largest_palindrome)  
