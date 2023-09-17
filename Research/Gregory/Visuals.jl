#=
# Collection of functions for visuallizing properties of integrators
#
# Author: Caleb Jacobs
# DLM: September 12, 2023
=#

include("GenGreg.jl")
using PlotlyJS
# using Plots
using Colors
using BenchmarkTools

"""
    plotPath(idx, g)
Display the path defined by idx through the grid g.
"""
function plotPath(idx::Vector{Int64}, g::Grid)
    x⃗ = real(g.z)
    y⃗ = imag(g.z)
    z⃗ = zeros(length(g.z))

    z⃗[idx] = [1 : length(idx)...]

    plt = plot(heatmap(
            x = x⃗, y = y⃗, z = z⃗,
            zsmooth = "none",
            colorscale = colors.viridis),
        Layout(width = 800, heights = 800))

    return plt
end

function plotPath(idx::Vector{Any}, g::Grid)
    x⃗ = real(g.z)
    y⃗ = imag(g.z)
    z⃗ = zeros(length(g.z))

    for path ∈ idx
        if typeof(path) <: Vector
            z⃗[path] = [1 : length(path)...]
        elseif typeof(path) <: Number
            z⃗[path] = 1
        end
    end

    plt = plot(heatmap(
            x = x⃗, y = y⃗, z = z⃗,
            zsmooth = "none",
            colorscale = colors.viridis),
        Layout(width = 800, heights = 800))

    return plt
end

"""
    plotGrid(g)

Plot grid `g` highlighting internal nodes, external nodes, and padding nodes.
"""
function plotGrid(g::Grid)
    x⃗ = real(g.z)
    y⃗ = imag(g.z)
    z⃗ = zeros(length(g.z))
    z⃗[g.e] .= 2
    z⃗[g.i] .= 1

    plt = plot(heatmap(
            x = x⃗, y = y⃗, z = z⃗,
            zsmooth = "none",
            zmin = 0, zmax = 2,
            colorscale = colors.viridis),
        Layout(width = 800, heights = 800))
end

function complexAbsPlot(z⃗, f⃗)
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    z⃗ = log10.(abs.(f⃗))
    z⃗[isinf.(z⃗)] .= -17

    display(findall(isnan.(z⃗)))
    
    plt = plot(heatmap(
            x = x⃗,
            y = y⃗,
            z = z⃗,
            zsmooth = "none",
            zmin = -16, zmax = 1,
            colorscale = colors.viridis),
        Layout(
            width = 800, height = 800))
    
    return plt
end

function complexPlot(z⃗, f⃗)
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    args = angle.(f⃗[:])
    args[args .< 0] .+= 2π

    plt = plot(heatmap(
            x = x⃗,
            y = y⃗,
            z = args,
            zsmooth = "best",
            zmin = 0, zmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["-π", "0", "π"])),
        Layout(
            width = 800, height = 800))
    
    return plt
end

function complexPlot3d(z⃗, f⃗)
    x⃗ = real(z⃗)
    y⃗ = imag(z⃗)
    z⃗ = abs.(f⃗)
    args = angle.(f⃗)
    args[args .< 0] .+= 2π

    plt = plot(surface(
            x = x⃗, y = y⃗, z = z⃗,
            surfacecolor = args,
            cmin = 0, cmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["-π", "0", "π"])))
        
    return plt
end
