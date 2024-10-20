function complex_phase_plot(z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    return heatmap(z.x, z.y, θ;
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   kwargs...
                  )
end

function complex_phase_plot!(ax::Union{Axis,Scene}, z::ComplexGrid, f::Matrix{<: Complex}; kwargs...)
    θ = mod.(angle.(f), 2π)

    return heatmap!(ax, z.x, z.y, θ;
                   colormap = :hsv,
                   colorrange = (0, 2π),
                   kwargs...
                  )
end
