include("sphericalharmonics.jl")

vand(nodes::Vector{<:Tuple}, deg) = [sphericalharmonic(x, n) for n in 0:(deg + 1)^2 - 1, x in nodes]
