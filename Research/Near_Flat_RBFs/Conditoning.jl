#=
# Experimenting with the condition number of collocation matrix
#
# Author: Caleb Jacobs
# DLM: 20-04-2022
=#

using Plots
using LinearAlgebra
using DoubleFloats

default(levels = 100, aspect_ratio = :equal, c = :viridis)

Φ(r, ε) = ℯ^(-(ε * r)^2)

"""
    colloc(X, ϕ, ε)

Construct the collocation matrix given nodes `X`, radial function `ϕ`, and shape paremeter `ε`.
"""
function colloc(X, ϕ, ε)
    N = size(X, 1)

    return [ϕ(norm(X[:,i] - X[:,j]), ε) for i ∈ 1:N, j ∈ 1:N]
end

"""
    getConds(X, ϕ, ε)

Compute the condition numbers of the collocation matrix at various values of the shape paremeter.
"""
function getConds(X, ϕ, ε)
    N = length(ε)

    conds = zeros(Float64, size(ε))

    for i ∈ 1:N
        conds[i] = cond(colloc(X, ϕ, ε[i]))
    end

    return conds
end


"""
    getε(rN, iN, xa, xb, ya, yb)

Generate complex grid of values of shape parameter `ε`.
"""
function getε(rN, iN, xa = 0, xb = 1, ya = 0, yb = 1)
    x = [range(-.5, .5, length = rN)...]
    y = [range(-.5, .5, length = iN)...]

    X = repeat([x...]', iN, 1)
    Y = repeat(y, 1, rN)

    ε = X + im * Y

    return ε
end

"""
    randDisk(N)

Generate `N` random nodes in a unit disk.
"""
function randDisk(N)
    Random.seed!(315)
    r = sqrt.(rand(N))
    r[1] = 0
    θ = 2π * rand(N)

    x = r .* cos.(θ)
    y = r .* sin.(θ)

    XY = [x'; y']
end

"""
    getNodes(N, κ)

Generate nodes over surface of a circle with curvature κ
"""
function getNodes(N, κ)
    x = randDisk(N)

    if κ > 0
        z = [sqrt(κ^(-2) - norm(x[:, i]).^2) - κ^(-1) for i ∈ 1:N]
    else
        z = zeros(N)
    end

    display(norm(z))

    X = [x; z']
end

"""
    genPlot(ϕ, N, rN, iN, κ)

Generate condition number plot of collocation with specified parameters.
"""
function genPlot(ϕ, N, rN, iN, κ)
    ε = getε(rN, iN)
    X = getNodes(N, κ)

    conds = getConds(X, ϕ, ε)

    plot(real.(ε[1,:]), imag.(ε[:,1]), log10.(conds), st = contour)
end

function streamPlots(ϕ, N, rN, iN, κ)
    for i = 1 : length(κ)
        display(genPlot(ϕ, N, rN, iN, κ[i]))
        sleep(.1)
    end
end

κ = 1 ./ (1.5 .^ [0:40...])
streamPlots(Φ, 101, 100, 100, κ)
