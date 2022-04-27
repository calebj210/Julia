#=
# Experimenting with the condition number of collocation matrix
#
# Author: Caleb Jacobs
# DLM: 26-04-2022
=#

using Plots
using Plots.PlotMeasures
using LinearAlgebra
using DoubleFloats
using Random
using LaTeXStrings
using Printf

default(levels = 300, aspect_ratio = :none, c = :viridis)

Φ(r, ε) = ℯ^(-(ε * r)^2)

"""
    colloc(X, ϕ, ε)

Construct the collocation matrix given nodes `X`, radial function `ϕ`, and shape paremeter `ε`.
"""
function colloc(X, ϕ, ε)
    N = size(X, 2)

    return [ϕ(norm(X[:,i] - X[:,j]), ε) for i ∈ 1:N, j ∈ 1:N]
end

"""
    hexGen(N; minx, maxx, miny, maxy)

Generate a hexagonal grid with approximately 'N' nodes that spans the rectangle
specified by 'minx', 'maxx', 'miny', and 'maxy'

See 'hexBand'.
"""
function hexGen(N; minx = -1, maxx = 1, miny = -1, maxy = 1)
    # Compute bound ranges
    Δx = maxx - minx
    Δy = maxy - miny
    
    # Compute number of nodes in each dimension
    n = round(Int,sqrt(N*Δx/Δy))
    m = round(Int, sqrt(N*Δy/Δx))
    
    # Compute base grid
    x = 0:n-1
    y = 0:m-1
    nodes = [repeat(x, inner = (m,1))' / ((n-1)+0.5) * Δx;
             repeat(y, outer = (n,1))' / (m-1)*Δy]

    # Compute x offsets
    xOff = repeat(repeat([0,0.5 / ((n-1)+0.5) * Δx],
                         outer = (ceil(Int,m/2),1))[1:m],
                  outer = (n,1)) .+ minx

    # Offset the nodes
    nodes[1,:] += xOff
    nodes[2,:] .+= miny
    
    return nodes
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
        # conds[i] = 1 / minimum(abs.(eigvals(colloc(X, ϕ, ε[i]))))
    end

    return conds
end


"""
    getε(rN, iN, xa, xb, ya, yb)

Generate complex grid of values of shape parameter `ε`.
"""
function getε(rN, iN, xa = 0, xb = 1, ya = 0, yb = 1)
    x = [range(0, 1, length = rN)...]
    y = [range(0, 1, length = iN)...]

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
    Random.seed!(210)
    r = sqrt.(rand(N))
    θ = 2π * rand(N)

    x = r .* cos.(θ)
    y = r .* sin.(θ)

    XY = [x'; y']
end

"""
    hexDisk(N)

Generate hexagonal node set of size less than `N` imposed in a unit disk.
"""
function hexDisk(N)
    xtmp = hexGen(1.5N) 
    N = size(xtmp, 2)

    xFull = zeros(2, N)
    j = 1
    for i = 1:N
        if norm(xtmp[:, i]) ≤ 1
            xFull[:, j] = xtmp[:, i]
            j += 1
        end
    end

    @printf "Number of nodes n = %d\n" N

    return xFull[:, 1 : j - 1]
end

"""
    getNodes(N, κ)

Generate nodes over surface of a circle with curvature κ
"""
function getNodes(N, κ)
    # x = randDisk(N)
    # N = size(x, 2)
    x = range(-1, 1, length = N)'

    if κ > 0
        z = [sqrt(κ^(-2) - norm(x[:,i]).^2) - κ^(-1) for i ∈ 1:N]
    else
        z = Double64.(zeros(N))
    end

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
    # scatter(real.(ε[:]), imag.(ε[:]), log10.(conds[:]), ms = 2, mα = 0.5, msα = 0)
end

function streamPlots(ϕ, N, rN, iN, κ)
    for i = 1 : length(κ)
        display(genPlot(ϕ, N, rN, iN, κ[i]))
    end
end

function test1()
    κ = 1 ./ (1.5 .^ [0:10...])
    streamPlots(Φ, 9, 75, 75, Double64.(κ))
end

function plotData(N)
    κ = [0, 0.5, 1]
    l = @layout [a b c]

    p = plot(aspect_ratio = :none, layout = l, size = (1000, 300))
    for i ∈ 1:3
        X = getNodes(N, κ[i])

        scatter!(X[1,:], X[2,:], X[3,:],
                 subplot = i, 
                 xlims = (-1.1, 1.1),
                 ylims = (-1.1, 1.1),
                 zlims = (-1.1, 1.1),
                 color = :black,
                 ms = 1,
                 legend = false,
                 dpi = 90,
                 title = L"\kappa = %$(κ[i])",
                 xlabel = "x",
                 ylabel = "y",
                 zlabel = "z",
                 aspect_ratio = :equal,
                 bottom_margin = 10mm)
    end

    return p
end

function discreteCond(ϕ, N)
    κ = Double64.([1, 0.75, 0.5, 0.25, 0.1, 0])
    l = @layout [a b c; d e f]

    p = plot(aspect_ratio = :none, layout = l, size = (1200, 500))
    for i ∈ 1 : length(κ)
        ε = getε(100, 100)
        X = getNodes(N, κ[i])
        
        conds = getConds(X, ϕ, ε)

        plot!(real.(ε[1,:]), imag.(ε[:,1]), log10.(conds), 
              st = contour, 
              subplot = i,
              title = L"\kappa = %$(κ[i])",
              colorbar_title = L"\log_{10}(\mathrm{cond}(A))",
              xlabel = L"Re(\varepsilon)",
              ylabel = L"Im(\varepsilon)",
              xlims  = (0,1),
              ylims  = (0,1),
              aspect_ratio = :equal,
              right_margin = 6mm,
              left_margin  = 6mm)
    end

    return p
end
