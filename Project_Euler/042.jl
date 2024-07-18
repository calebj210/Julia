words = String.(strip.(split(read("InputFiles/042_words.txt", String), ','), '"'))

function wordvalue(name::String)
    return sum(Int.(collect(name))) - length(name) * (Int('A') - 1)
end
    
function istriangular(n::Int)
    return isinteger((sqrt(8n + 1) - 1) / 2)
end
istriangular(s::String) = istriangular(wordvalue(s))

function count_triangular_words()
    count(istriangular.(words))
end

count_triangular_words()
