#=
# Plotting routines for FD stencils
#
# Author: Caleb Jacobs
# DLM: April 10, 2025
=#
using CairoMakie
using LaTeXStrings
using LinearAlgebra

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

    fig = Figure()
    ax = Axis(fig[1,1], 
        xlabel = "Distance from Origin",
        ylabel = L"|\omega|",
        yscale = log10
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

    Legend(fig[1,2], ax)

    resize_to_layout!(fig)

    return fig
end
