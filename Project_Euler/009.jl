function find_pythagorean_triplet_with_sum(sum_value)
    for a in 1:sum_value ÷ 3 # a cannot be more than 1/3 of the sum
        b = (sum_value^2 - 2 * a * sum_value) ÷ (2 * sum_value - 2 * a)
        c = sum_value - a - b
        if a^2 + b^2 == c^2
            return prod([a, b, c])
        end
    end
    return nothing  # No triplet found
end
