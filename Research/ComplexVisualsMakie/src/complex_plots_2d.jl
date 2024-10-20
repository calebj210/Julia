function complex_phase_plot(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    return heatmap(z.x, z.y, θ;
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   kwargs...
                  )
end

function complex_phase_plot!(ax::Axis, z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    return heatmap!(ax, z.x, z.y, θ;
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   kwargs...
                  )
end

function complex_phase_plot!(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    return heatmap!(z.x, z.y, θ;
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   kwargs...
                  )
end

function complex_color_wheel()
    fig = Figure(backgroundcolor = :transparent)
    ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)

    z = complex_grid(1, 500)

    f = collect(z)
    f[abs.(z) .>= 1] .= NaN

    plt = complex_phase_plot!(ax, z, f)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)

    return fig
end
