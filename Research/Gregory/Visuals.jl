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
    plt = plot(g.z, st = :scatter, c = :red, ratio = 1, legend = false,
          mα = 0.9, msw = 0)
    
    plot!(g.z[idx], st = :scatter, mz = 1 : length(idx), c = :viridis,
          mα = 0.9, msw = 0)

    return plt
end
function plotPath(idx::Tuple{Vector, Vector}, g::Grid)
    plt = plot(g.z, st = :scatter, c = :red, ratio = 1, legend = false,
               mα = 0.9, msw = 0)
    
    for i ∈ idx
        if !isempty(i)
            plot!(g.z[i], st = :scatter, mz = 1 : length(i), c = :viridis,
                  mα = 0.9, msw = 0)
        end
    end

    return plt
end

"""
    plotGrid(g)

Plot grid `g` highlighting internal nodes, external nodes, and padding nodes.
"""
function plotGrid(g::Grid)
    plt = plot(g.z[g.i], st = :scatter, c = :green, mα = 0.75, msw = 0, ms = 4.5, label = "Taylor")
    plot!(g.z[g.e], st = :scatter, c = :blue,       mα = 0.75, msw = 0, ms = 4.5, label = "Gregory")
    plot!(g.z[g.p], st = :scatter, c = :red,        mα = 0.75, msw = 0, ms = 4.5, label = "Padding")
    plot!(ratio = 1, dpi = 300)
end

# function plotComplexErrors(tru, Df, g::Grid)
#     plt = plot(g.z, st = :scatter, 
#                c   = :viridis,
#                mz  = log10.(abs.(tru - Df)),
#                mα  = 0.75,
#                msw = 0,
#                ratio = 1,
#                clims = (-16, 1),
#                label = false)
# 
#     return plt
# end

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
