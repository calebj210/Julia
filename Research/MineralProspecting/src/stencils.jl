using CairoMakie, LaTeXStrings

function stencil_convergence_2D!(ax)
    hs = 10 .^ (0:-0.5:-7)
    FD =  [
        -1 0 1;
        -1 0 1;
        -1 0 1
    ] / 6
    HFD = [
        -1 0 1;
        -4 0 4;
        -1 0 1
    ] / 12

    F((x,y)) = (x + 2) / ((x + 2)^2 + y^2)
    dF = -1/4

    err_fd = Vector{Float64}(undef, length(hs))
    err_hfd = Vector{Float64}(undef, length(hs))
    for (idx, h) in enumerate(hs)
        xs = [(x*h,y*h) for y in -1:1, x in -1:1]
        f = F.(xs)
        err_fd[idx] = abs(sum(f .* FD) / h - dF) / abs(dF)
        err_hfd[idx] = abs(sum(f .* HFD) / h - dF) / abs(dF)
    end

    lines!(ax, hs, err_fd, label = "FD")
    lines!(ax, hs, err_hfd, label = "HFD")

    lines!(ax, hs, hs.^2, linestyle = :dash, label = L"\mathcal{O}(h^2)")
    lines!(ax, hs, hs.^4, linestyle = :dash, label = L"\mathcal{O}(h^4)")

    axislegend(ax, position = :rt)
end

function stencil_convergence_3D!(ax)
    hs = 10 .^ (0:-0.5:-7)
    FD = zeros(Float64, (3,3,3))
    FD[:,:,1] = [
        1 1 1;
        1 1 1;
        1 1 1
    ] / -18
    FD[:,:,2] = [
        0 0 0;
        0 0 0;
        0 0 0
    ]
    FD[:,:,3] = [
        1 1 1;
        1 1 1;
        1 1 1
    ] / 18
    
    HFD = zeros(Float64, (3,3,3))
    HFD[:,:,1] = [
        0 1 0;
        1 2 1;
        0 1 0
    ] / -12
    HFD[:,:,2] = [
        0 0 0;
        0 0 0;
        0 0 0
    ]
    HFD[:,:,3] = [
        0 1 0;
        1 2 1;
        0 1 0
    ] / 12

    F((x,y,z)) = (z + 2) / sqrt(x^2 + y^2 + (z + 2)^2)^3
    dF = -1/4

    err_fd = Vector{Float64}(undef, length(hs))
    err_hfd = Vector{Float64}(undef, length(hs))
    for (idx, h) in enumerate(hs)
        xs = [(x*h,y*h,z*h) for y in -1:1, x in -1:1, z in -1:1]
        f = F.(xs)
        err_fd[idx] = abs(sum(f .* FD) / h - dF) / abs(dF)
        err_hfd[idx] = abs(sum(f .* HFD) / h - dF) / abs(dF)
    end

    lines!(ax, hs, err_fd, label = "FD")
    lines!(ax, hs, err_hfd, label = "HFD")

    lines!(ax, hs, hs.^2, linestyle = :dash, label = L"\mathcal{O}(h^2)")
    lines!(ax, hs, hs.^4, linestyle = :dash, label = L"\mathcal{O}(h^4)")

    axislegend(ax, position = :rt)
end

function stencil_convergence()
    set_theme!(theme_latexfonts())
    fig = Figure()
    ax2d = Axis(
        fig[1,1],
        title = L"2D $\partial / \partial y$",
        xlabel = L"Grid Spacing $h$",
        ylabel = "Relative Error",
        xreversed = true,
        xscale = log10,
        yscale = log10,
        limits = (nothing, (1e-16, 1e1)),
    )
    ax3d = Axis(
        fig[1,2],
        title = L"3D $\partial / \partial z$",
        xlabel = L"Grid Spacing $h$",
        ylabel = "Relative Error",
        xreversed = true,
        xscale = log10,
        yscale = log10,
        limits = (nothing, (1e-16, 1e1)),
    )

    stencil_convergence_2D!(ax2d)
    stencil_convergence_3D!(ax3d)

    rowsize!(fig.layout, 1, Aspect(2,1))

    resize_to_layout!(fig)

    return fig
end
