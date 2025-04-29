#= 
# Code to show the interaction of various 2F1 transformations
#
# Author: Caleb Jacobs
# DLM: April 28, 2025
=#

using CairoMakie, ComplexVisuals
include("pFq.jl")

function failed_mapping(tru = nothing)
    a, b, c = (-19.139483726210333, -14.134535062801518, -24.30112282763963)
    z = complex_square_grid(2, 300)

    if isnothing(tru)
        tru = johansson_2f1.(a, b, c, z)
    end
    f = taylor_2f1.(a, b, c, z)

    errs = clean_error.(f,tru)

    fig = Figure()

    ax0 = Axis(fig[1,1],
        title = L"z \mapsto z",
        limits = ((-3,3), (-3,3)),
    )
    ax1 = Axis(fig[1,2],
        title = L"z \mapsto \frac{1}{z}",
        limits = ((-3,3), (-3,3)),
    )
    ax2 = Axis(fig[1,3],
        title = L"z \mapsto 1 - z",
        limits = ((-3,3), (-3,3)),
    )
    ax3 = Axis(fig[2,1],
        title = L"z \mapsto 1 - \frac{1}{z}",
        limits = ((-3,3), (-3,3)),
    )
    ax4 = Axis(fig[2,2],
        title = L"z \mapsto \frac{1}{1 - z}",
        limits = ((-3,3), (-3,3)),
    )
    ax5 = Axis(fig[2,3],
        title = L"z \mapsto \frac{z}{z - 1}",
        limits = ((-3,3), (-3,3)),
    )

    contour!(ax0, reim(z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax1, reim(1 ./ z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax2, reim(1 .- z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax3, reim(1 .- 1 ./ z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax4, reim(1 ./ (1 .- z))..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax5, reim(z ./ (z .- 1))..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )

    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 2, Aspect(1, 1))
    colsize!(fig.layout, 3, Aspect(1, 1))

    # Colorbar(fig[1,3], ax0)

    resize_to_layout!(fig)

    return (fig, tru)
end

function mapped_regions()
    z = -2 + 2im .+ cispi.(range(0, 2, 100))

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax0 = Axis(fig[1,1],
        title = L"z \mapsto z",
        limits = ((-3,3), (-3,3)),
    )
    ax1 = Axis(fig[1,2],
        title = L"z \mapsto \frac{1}{z}",
        limits = ((-3,3), (-3,3)),
    )
    ax2 = Axis(fig[1,3],
        title = L"z \mapsto 1 - z",
        limits = ((-3,3), (-3,3)),
    )
    ax3 = Axis(fig[2,1],
        title = L"z \mapsto 1 - \frac{1}{z}",
        limits = ((-3,3), (-3,3)),
    )
    ax4 = Axis(fig[2,2],
        title = L"z \mapsto \frac{1}{1 - z}",
        limits = ((-3,3), (-3,3)),
    )
    ax5 = Axis(fig[2,3],
        title = L"z \mapsto \frac{z}{z - 1}",
        limits = ((-3,3), (-3,3)),
    )

    lines!(ax0, reim(z)..., 
        color = 1:length(z)
    )
    lines!(ax1, reim(1 ./ z)...,
        color = 1:length(z),
    )
    lines!(ax2, reim(1 .- z)...,
        color = 1:length(z),
    )
    lines!(ax3, reim(1 .- 1 ./ z)...,
        color = 1:length(z),
    )
    lines!(ax4, reim(1 ./ (1 .- z))...,
        color = 1:length(z),
    )
    lines!(ax5, reim(z ./ (z .- 1))...,
        color = 1:length(z),
    )
    
    return fig
end

function clean_error(f,t)
    err = abs.((f - t) ./ t)

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-16
        err = 1e-16
    end

    return err
end
