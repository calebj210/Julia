names = String.(strip.(split(read("InputFiles/022_names.txt", String), ','), '"'))

function nameScore(name::String)
    return sum(Int.(collect(name))) - length(name) * (Int('A') - 1)
end

function nameScoreSum(list::Vector{String})
    list = sort(list)

    s = 0
    for (n, name) ∈ enumerate(list)
        s += nameScore(name) * n
    end

    return s
end
