#=
# Collection of functions for visuallizing properties of integrators
#
# Author: Caleb Jacobs
# DLM: November 1, 2023
=#

include("GenGreg.jl")
using PlotlyJS
using Colors

"""
    plotPath(idx, g)
Display the path defined by idx through the grid g.
"""
function plotPath(idx::Vector{Path}, g::Grid)
    x⃗ = real(g.z)
    y⃗ = imag(g.z)
    z⃗ = zeros(length(g.z))

    for p ∈ idx
        z⃗[p.i] = 1
        z⃗[p.p] = [1 : length(p.p)...]
        z⃗[p.f] = 1
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
    z⃗[g.e] .= 1
    z⃗[g.i] .= 2
    z⃗[g.ib] = 2 * [0 : (length(g.ib) - 1)...] / length(g.ib)

    plt = plot(heatmap(
            x = x⃗, y = y⃗, z = z⃗,
            zsmooth = "none",
            zmin = 0, zmax = 2,
            colorscale = colors.viridis),
        Layout(width = 800, heights = 800))
end

function complexAbsPlot(z⃗, f⃗; logscale = false, title = "")
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    if logscale
        z⃗ = log10.(abs.(f⃗))
        z⃗[isinf.(z⃗)] .= -17
    else
        z⃗ = abs.(f⃗)
    end

    layout = Layout(
        width = 800, height = 800,
        xaxis = attr(
            title = attr(
                text = "Re(z)",
                font_size = 20
            )
        ),
        yaxis = attr(
            title = attr(
                text = "Im(z)",
                font_size = 20
            )
        ),
        title = attr(
            text = title,
            font_size = 25,
            y = 0.96,
            x = 0.5,
            xanchor = "center",
            yanchor = "top"
        )
    )

    if logscale
        plt = plot(heatmap(
                x = x⃗,
                y = y⃗,
                z = z⃗,
                zsmooth = "none",
                zmin = -16, zmax = 1,
                colorscale = colors.viridis),
               layout)
    else
        plt = plot(heatmap(
                x = x⃗,
                y = y⃗,
                z = z⃗,
                zsmooth = "none",
                colorscale = colors.viridis),
               layout)
    end
    
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
            zsmooth = "none",
            zmin = 0, zmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["0", "π", "2π"])),
            Layout(width = 800, height = 800))
    
    return plt
end

function complexPlot3d(z⃗::Matrix, f⃗::Matrix...; T = 1)
    x⃗ = real(z⃗)
    y⃗ = imag(z⃗)
        
    plts = Vector{GenericTrace}()

    for p ∈ f⃗
        if T == 2
            z⃗ = real(p)
            args = angle.(p)
            args[args .< 0] .+= 2π
        elseif T == 3
            z⃗ = imag(p)
            args = angle.(p)
            args[args .< 0] .+= 2π
        else
            z⃗ = abs.(p)
            args = angle.(p)
            args[args .< 0] .+= 2π
        end

        push!(plts, surface(
                        x = x⃗, y = y⃗, z = z⃗,
                        surfacecolor = args,
                        cmin = 0, cmax = 2π,
                        colorscale = colors.hsv,
                        colorbar = attr(
                            tickmode = "array",
                            tickvals = [0, π, 2π],
                            ticktext = ["-π", "0", "π"])))
    end

    if T == 2
        title = "Re(f)"
    elseif T == 3
        title = "Im(f)"
    else
        title = "Abs-Arg(f)"
    end

    layout = Layout(
        width = 800, height = 800,
        title = attr(
            text = title,
            font_size = 25,
            y = 0.96,
            x = 0.5,
            xanchor = "center",
            yanchor = "top"
        )
    )

    plt = plot(plts, layout)

    return plt
end

function complexPlot3d(z::Vector, f::Vector...; T = 1)
    Z = reshape(z, round(Int64, sqrt(length(z))), :)

    N = round(Int64, sqrt(length(f[1])))

    F = reshape.(f, N , :)

    return complexPlot3d(Z, F..., T = T)
end
