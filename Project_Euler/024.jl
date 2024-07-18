using Combinatorics

function getnthperm(n)
    return parse(Int, join(nthperm(0:9, n)))
end

getnthperm(10^6)
