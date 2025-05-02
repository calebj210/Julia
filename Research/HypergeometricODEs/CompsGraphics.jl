using CairoMakie, ComplexVisuals
include("pFq.jl")

pslog(z) = sign(z) * log10(abs(z) + 1)

function Levin_Test_Plot()
    a,b,c = (1,-9/2,-9/4)
    z = complex_square_grid(10, 300)
    f = taylor_2f1.(a, b, c, z)

    set_theme!(theme_latexfonts())
    tickvals = 
    ticks = (pslog.([0,1e2,1e4,1e6]), [L"0", L"10^2", L"10^4", L"10^6"])
    fig = Figure()
    ax = Axis3(fig[1,1],
        title = L"Abs-Arg of ${_2}F_1(1, -\frac{9}{2}; -\frac{9}{4}; z)$",
        titlesize = 20,
        xlabelsize = 15,
        ylabelsize = 15,
        zlabelsize = 15,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"pseudolog$_{10}(f)$",
        zticks = ticks,
        elevation = 0.2pi,
    )
    complexsurface!(ax, z, pslog.(f))

    complex_color_wheel!(ax.scene, radius = .5, center = (8,3.5))

    return fig
end
