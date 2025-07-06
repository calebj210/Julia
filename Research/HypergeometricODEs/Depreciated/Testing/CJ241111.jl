#=
# Timing tests for my CJ241111 handout
#
# Author: Caleb Jacobs
# DLM: November 11, 2024
=#

using ComplexVisuals, CairoMakie
using BenchmarkTools
include("../pFq.jl")

function time_test()
    z = complex_square_grid(10, 300)

    print("Running fixed step size...")
    fixed = average_time.(1, -9/2, -9/4, z; H = 0.1)
    println(" done.")

    print("Running optimal step size...")
    adapt = average_time.(1, -9/2, -9/4, z)
    println(" done.")

    set_theme!(merge(complex_theme, theme_latexfonts()))
    fig = Figure()
    ax1 = Axis(fig[1,1], title = "Fixed Step Size Time", xlabel = "")
    ax2 = Axis(fig[1,2], title = "Optimal Step Size Time", xlabel = "", ylabel = "")
    ax3 = Axis(fig[2,1], title = "Fixed Step Size Error")
    ax4 = Axis(fig[2,2], title = "Optimal Step Size Error", ylabel = "")

    tmin = min(minimum(fixed), minimum(adapt))
    tmax = min(maximum(fixed), maximum(adapt))
    trange = (tmin, tmax)

    plt1 = heatmap!(ax1, z, fixed; colorscale = log10, colorrange = trange)
    plt2 = heatmap!(ax2, z, adapt; colorscale = log10, colorrange = trange)
    Colorbar(fig[1,3], plt1)

    tru = johansson_2f1.(1, -9/2, -9/4, z, bits = 256)

    fixed = taylor_2f1.(1, -9/2, -9/4, z; H = 0.1)
    fixed = abs.((fixed - tru) ./ tru)

    adapt = taylor_2f1.(1, -9/2, -9/4, z)
    adapt = abs.((adapt - tru) ./ tru)

    emin = min(minimum(fixed), minimum(adapt))
    emin = emin == 0 ? 1e-17 : emin
    emax = min(maximum(fixed), maximum(adapt))
    erange = (emin, emax)
    display(erange)

    plt3 = heatmap!(ax3, z, fixed; colorscale = log10, colorrange = erange)
    plt4 = heatmap!(ax4, z, adapt; colorscale = log10, colorrange = erange)
    Colorbar(fig[2,3], plt3)

    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 2, Aspect(1, 1))

    resize_to_layout!(fig)

    return fig
end

average_time(a, b, c, z; N = 5, H = Inf) = median([@elapsed taylor_2f1(a, b, c, z; H = H) for i ∈ 1:N])
