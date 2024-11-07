function complex_color_wheel()
    fig = Figure(backgroundcolor = :transparent)
    ax = Axis(fig[1, 1], backgroundcolor = :transparent)
    hidespines!(ax)
    hidedecorations!(ax)

    z = complex_square_grid(1, 500)

    f = collect(z)
    f[abs.(z) .>= 1] .= NaN

    plt = phase!(ax, z, f)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)

    return fig
end

function complex_color_wheel!(scene::Scene; center = (.85,.15), radius = .1)
    relative = Makie.camrelative(scene)
    z = complex_square_grid(1, 500)

    f = collect(z)
    f[abs.(z) .>= 1] .= NaN

    plt = phase!(relative, z, f)
    scale!(plt, radius, radius)
    translate!(plt, center...)

    return plt
end
