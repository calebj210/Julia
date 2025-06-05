#= 
# Code to show the interaction of various 2F1 transformations
#
# Author: Caleb Jacobs
# DLM: May 27, 2025
=#

using CairoMakie, ComplexVisuals, LaTeXStrings
include("pFq.jl")

function sign_test(asign = -1, bsign = -1, csign = -1; tru = nothing, transformation = "")
    a, b, c = (asign * 19.139483726210333, bsign * 14.134535062801518, csign * 24.30112282763963)
    z = complex_square_grid(2, 300)

    if isnothing(tru)
        tru = johansson_2f1.(a, b, c, z)
    end
    f = taylor_2f1.(a, b, c, z)

    errs = clean_error.(f, tru)

    ticks = ( -16:4:0, [latexstring("10^{", p, "}") for p ∈ -16:4:0])
    set_theme!(theme_latexfonts())
    fig1 = Figure()
    fig2 = Figure()
    ax1 = Axis3(fig1[1,1],
        title = latexstring(
            "\${_2}F_1(", 
                (asign == -1 ? "-a," : "a,"), 
                (bsign == -1 ? "-b;" : "b;"), 
                (csign == -1 ? "-c;" : "c;"), 
            "z)\$ via ",
            transformation
        ),
        titlesize = 20,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        limits = (nothing, nothing, (-17,2)),
        zticks = ticks,
        protrusions = (60, 30, 30, 30)
    )
    ax2 = Axis3(fig2[1,1],
        title = L"Abs-Arg of ${_2}F_1(a,b;c;z)$",
        titlesize = 20,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"pseudolog$_{10}|f(z)|$",
        zlabeloffset = 30,
    )
    p1 = surface!(ax1, reim(z)..., log10.(errs), colorrange = (-17,1))
    complexsurface!(ax2, z, Makie.pseudolog10.(tru))

    Colorbar(fig1[1,2], p1, ticks = ticks)
    resize_to_layout!(fig1)

    return (fig1, fig2, tru)
end

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
        ylabel = L"\mathrm{Re}(z)",
        limits = ((-3,3), (-3,3)),
    )
    ax1 = Axis(fig[1,2],
        title = L"z \mapsto \frac{1}{z}",
        limits = ((-3,3), (-3,3)),
    )
    ax2 = Axis(fig[2,1],
        title = L"z \mapsto \frac{1}{1 - z}",
        ylabel = L"\mathrm{Re}(z)",
        limits = ((-3,3), (-3,3)),
    )
    ax3 = Axis(fig[2,2],
        title = L"z \mapsto 1 - z",
        limits = ((-3,3), (-3,3)),
    )
    ax4 = Axis(fig[3,1],
        title = L"z \mapsto 1 - \frac{1}{z}",
        xlabel = L"\mathrm{Im}(z)",
        ylabel = L"\mathrm{Re}(z)",
        limits = ((-3,3), (-3,3)),
    )
    ax5 = Axis(fig[3,2],
        title = L"z \mapsto \frac{z}{z - 1}",
        xlabel = L"\mathrm{Im}(z)",
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
    contour!(ax2, reim(1 ./ (1 .- z))..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax3, reim(1 .- z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax4, reim(1 .- 1 ./ z)..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )
    contour!(ax5, reim(z ./ (z .- 1))..., errs,
        colorscale = log10,
        colorrange = (1e-16, 1e1),
        levels = 10.0 .^ (-16:1),
    )

    # colsize!(fig.layout, 1, Aspect(1, 1))
    # colsize!(fig.layout, 2, Aspect(1, 1))
    # colsize!(fig.layout, 3, Aspect(1, 1))

    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Aspect(1, 1))
    rowsize!(fig.layout, 3, Aspect(1, 1))

    ticks = (-15:3:0, [latexstring("10^{", p, "}") for p ∈ -15:3:0])
    Colorbar(fig[1,3], limits = (-16, 1), colormap = :viridis, ticks = ticks, label = "Relative Error", labelrotation = -π/2)
    Colorbar(fig[2,3], limits = (-16, 1), colormap = :viridis, ticks = ticks, label = "Relative Error", labelrotation = -π/2)
    Colorbar(fig[3,3], limits = (-16, 1), colormap = :viridis, ticks = ticks, label = "Relative Error", labelrotation = -π/2)

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
