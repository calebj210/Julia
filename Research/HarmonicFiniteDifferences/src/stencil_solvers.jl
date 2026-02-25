include("harmonics.jl")

vand(nodes::Vector{<:NTuple{3,<:Number}}; deg = 0, N = (deg + 1)^2 - 1) =
    [sphericalharmonic(x, n) for n in 0:N, x in nodes]
vand(nodes::Vector{<:NTuple{2,<:Number}}; deg = 0, N = 2deg) =
    [circularharmonic(x, n) for n in 0:N, x in nodes]

function FD_weights(nodes; deg = 0, deriv = (3,1))
    A = vand(nodes; deg)
    b = zeros(BigFloat, size(A)[1])
    b[deriv[1]] = deriv[2]

    w = A \ b

    return w
end
