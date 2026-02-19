include("harmonics.jl")

vand(nodes::Vector{<:NTuple{3,<:Number}}; deg = 0, N = (deg + 1)^2 - 1) =
    [sphericalharmonic(x, n) for n in 0:N, x in nodes]
vand(nodes::Vector{<:NTuple{2,<:Number}}; deg = 0, N = 2deg) =
    [circularharmonic(x, n) for n in 0:N, x in nodes]
