function spiraldiag(N::Integer)
    if N % 2 == 0
        error("Only works for odd sized matrices.")
    end

    diags = [1]

    for n ∈ 1 : N ÷ 2
        step = 2n
        push!(diags, (diags[end] .+ step * [1:4...])...)
    end

    return diags
end

sum(spiraldiag(1001))
