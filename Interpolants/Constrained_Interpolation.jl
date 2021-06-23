#=
# Gradient constrained optimization
# Author: Caleb Jacobs
# Date last modified: 04/01/2021
=#

using Plots

  ϕ(x1, x2, m) = abs(x1 - x2)^m
ϕ_x(x1, x2, m) = m*(x1 - x2) * abs(x1 - x2)^(m-2)

function main(m)
    nodes = [0, 3, 1, 2]
    N = 2

    A1  = [  ϕ(nodes[i], nodes[j], m) for i ∈ 1:N, j ∈ 1:2N]
    Ax  = [ϕ_x(nodes[i], nodes[j], m) for i ∈ 1:N, j ∈ 1:2N]
    A   = [A1; Ax]

    b1  = zeros(N)
    bx  = [-1.0, 1.0]
    b   = [b1; bx]

    λ = [A; Ax] \ [b; bx]

    display(λ)

    # Begin plotting result
    t = 0 : 0.01 : 3

    f = [sum([λ[j] * ϕ(t[i], nodes[j], m) for j ∈ 1:2N]) for i ∈ 1:lastindex(t)]

    plotA = plot(t, f,
                 ratio = 1)
    display(plotA)
end

main(3)
