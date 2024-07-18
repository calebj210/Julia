using Combinatorics

function find_pandigital_products()
    pandigital_products = Vector{Int}()

    for p ∈ permutations(1 : 9)
        prd = parse(Int, join(p[end - 3:end]))
        if prd ∈ pandigital_products
            continue
        end

        a = parse(Int, join(p[1]))
        b = parse(Int, join(p[2:5]))
        if a * b == prd
            push!(pandigital_products, prd)
            continue
        end

        a = parse(Int, join(p[1:2]))
        b = parse(Int, join(p[3:5]))
        if a * b == prd
            push!(pandigital_products, prd)
            continue
        end
    end

    return pandigital_products
end

sum(find_pandigital_products())
