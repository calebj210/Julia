isfactorialsum(n::Int) = sum(factorial.(digits(n))) == n

function get_factorial_sums(N::Int)
    factorial_sums = Vector{Int}()
    
    for n ∈ 3:N
        if isfactorialsum(n)
            push!(factorial_sums, n)
        end
    end

    return factorial_sums
end

sum(get_factorial_sums(10^5))
