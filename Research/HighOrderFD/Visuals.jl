#=
# Plotting routines for FD stencils
#
# Author: Caleb Jacobs
# DLM: May 2, 2025
=#
using CairoMakie
using LaTeXStrings
using LinearAlgebra

include("2DHarmonicFD.jl")

"""
    distancePlot2D(Z, D)
Create a log plot of the absolute value of the stencil weights in `D` as a function of distance from the origin for nodes in `Z`.
"""
function distancePlot2D(Z, D)
    ω = vec(Float64.(abs.(D)))
    z = vec(norm.(Z))

    ω[ω .== 0] .= ω[1]

    x = range(0, maximum(z), length = 100)
    y = maximum(ω) * exp.(-π / 2 * x.^2)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(fig[1,1], 
        xlabel = "Distance from Origin",
        ylabel = "|Weight|",
        # xlabelsize = 15,
        # ylabelsize = 15,
        yscale = log10,
    )

    lines!(ax, x, y,
        linestyle = :dash,
        color = :black,
        label = L"C e^{-\frac{x^2}{2}}",
    )

    scatter!(ax, z, ω,
        color = :gray,
        label = "Weights",
    )

    axislegend(ax, labelsize = 22)

    resize_to_layout!(fig)

    return fig
end


function Δdistanceplot()
    setprecision(1024)

    Z = get_centered_stencil(15)
    D = getΔWeights(15, 139)

    fig = distancePlot2D(Z, D)
    fig.content[1].title = L"$\Delta$ Stencil Weights"
    fig.content[1].titlesize = 24
    fig.content[1].xlabelsize = 22
    fig.content[1].ylabelsize = 22
    fig.content[1].xticklabelsize = 22
    fig.content[1].yticklabelsize = 22
    fig.content[1].xticks = (0:2:10)
    fig.content[1].yticks = (
        [1e-64, 1e-48, 1e-32, 1e-16, 1e0],
        [latexstring("10^{", p, "}") for p ∈ [-64, -48, -32, -16, 0]]
    )

    return fig
end
