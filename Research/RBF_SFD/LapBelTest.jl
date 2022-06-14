include("RBF_SFD.jl")

function lapTest(N, n, m, o)
    z(s) = s[1]^2 * s[2]^3
    f(s) = s[1]^2 + s[2]^2

    x = range(0, 1, length = N)
    nodes = [repeat(x, inner = length(x))';
             repeat(x, outer = length(x))';
             zeros(1, length(x)^2)]
    nodes[3, :] = [z(nodes[1:2, i]) for i ∈ 1:length(x)^2]

    F = [f(nodes[1:2, i]) for i ∈ 1:length(x)^2]

    coms = getCommons(nodes, n, m, o)
    # DΔ   = discΔ(nodes, coms) 

    # approx = DΔ * [f(nodes[1:2, i]) for i ∈ 1:length(x)^2]
    approx = computeΔ(nodes, F, coms)
    trues = [Δ(nodes[1:2, i], z, f) for i ∈ 1:length(x)^2]

    display(approx - trues)
    display(computeSurfDiv(nodes, F, coms, 1, 1))
end
