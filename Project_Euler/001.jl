function mult35(N::Integer)
    s = 0
    for i ∈ 1 : N - 1
        if i % 3 == 0 || i % 5 == 0
            s += i
        end
    end

    return s
end
