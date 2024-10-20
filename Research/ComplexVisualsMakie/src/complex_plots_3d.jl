function complex_surface_plot(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    fig = Figure()
    ax = Axis3(fig[1,1])
    plt = surface!(ax, z.x, z.y, abs.(f);
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   color = θ,
                   interpolate = false,
                   kwargs...
                  )
    ax.xlabel = "Re(z)"
    ax.ylabel = "Im(z)"
    ax.zlabel = "Abs(f)"

    return Makie.FigureAxisPlot(fig, ax, plt)
end

function complex_real_plot(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1,1])
    plt = wireframe!(ax, z.x, z.y, real.(f);
                   kwargs...
                  )
    ax.xlabel = "Re(z)"
    ax.ylabel = "Im(z)"
    ax.zlabel = "Re(f)"

    return Makie.FigureAxisPlot(fig, ax, plt)
end

function complex_imag_plot(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1,1])
    plt = wireframe!(ax, z.x, z.y, imag.(f);
                   kwargs...
                  )
    ax.xlabel = "Re(z)"
    ax.ylabel = "Im(z)"
    ax.zlabel = "Im(f)"

    return Makie.FigureAxisPlot(fig, ax, plt)
end
