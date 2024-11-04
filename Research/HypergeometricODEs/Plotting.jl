using Makie, ComplexVisualsMakie
import CairoMakie: heatmap

function heatmap(z::ComplexGrid, f::Matrix; kwargs...)
    fig, ax, plt = heatmap(z.x, z.y, f; kwargs...)

    ax.xlabel = "Re(z)"
    ax.ylabel = "Im(z)"

    return Makie.FigureAxisPlot(fig, ax, plt)
end

function grid_error_plot(z, f, tru; kwargs...)
    fig, ax, plt = heatmap(z, abs.((f - tru) ./ tru); 
                           colorscale = log10,
                           colorrange = (1e-16, 1e0),
                           kwargs...)

    colsize!(fig.layout, 1, Aspect(1,1))
    Colorbar(fig[1,2], plt)
    resize_to_layout!(fig)

    return Makie.FigureAxisPlot(fig, ax, plt)
end
