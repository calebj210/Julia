using CairoMakie, LaTeXStrings
using SparseArrays

function sparse_structure(x::Vector, deg)
    A = vand(x; deg)
    B = qr(A')
    D = diagm([1/norm((B.R')[n,:]) for n ∈ 1:size(B.R')[1]])
    idx = D * abs.(B.R') .> 1e-13
    C = spzeros(size(B.R'))
    C[idx] = log10.(abs.(B.R'[idx]))

    fig = Figure()
    ax = Axis(fig[1,1],
        yreversed = true,
        aspect = DataAspect(),
        xlabel = "nz = $(length(C.nzval))",
        title = "L",
    )
    plt = spy!(ax, C')

    Colorbar(fig[1,2], plt, label = L"\log_{10}|\text{val}|")

    colsize!(fig.layout, 1, Aspect(1,.9))
    resize_to_layout!(fig)

    return fig
end

function onesided_2d(;sizes = [(0,29), (1,9), (7,1)], deg = 14)
    fig = Figure()
    
    for (i,sz) in enumerate(sizes)
        nodes = onesided_nodes(sz..., type = BigFloat)

        w = FD_weights(nodes; deg, deriv = (3,1))

        ax = Axis3(
            fig[1,i],
            aspect = (15,30,30),
            limits = ((-7.5,7.5), (-0.5,30.5), (0,80)),
            azimuth = 0.6pi,
            elevation = .3,
            title = "$(2sz[1] + 1)x$(sz[2] + 1)",
            zlabel = L"|w|",
        )

        barplot3d!(
            ax, first.(nodes), last.(nodes), abs.(w), 
            colorrange = (1e-2,80),
            colormap = :viridis,
        )
    end

    Colorbar(fig[1,4], colorrange = (1e-2,80), colormap = :viridis, label = L"|w|")
    rowsize!(fig.layout, 1, Aspect(1,2))
    resize_to_layout!(fig)

    return fig
end

function plot_weights(nodes, w; axis = (;), barplot3d = (;))
    xrng = round(Int64, abs(-(extrema(first.(nodes))...))) + 1
    yrng = round(Int64, abs(-(extrema( last.(nodes))...))) + 1

    fig = Figure()
    
    ax = Axis3(
        fig[1,1];
        aspect = (xrng, yrng, max(xrng, yrng)),
        axis...
    )

    plt = barplot3d!(
        ax, first.(nodes), last.(nodes), w;
        colormap = :viridis,
        barplot3d...
    )

    Colorbar(fig[1,4], plt, label = L"|w|")
    resize_to_layout!(fig)

    return fig

end

function barplot3d(x, y, z; kwargs...)
    pos = collect(zip(x, y))
    sizes = [Vec3d(1, 1, z) for z in z]

    figaxplt = meshscatter(
        pos;
        markersize = sizes,
        marker = Rect3d((-0.5, -0.5, eps()), (1, 1, 1)),
        color = z,
        kwargs...
    )

    return figaxplt
end

function barplot3d!(ax, x, y, z; kwargs...)
    pos = collect(zip(x, y))
    sizes = [Vec3(1, 1, z) for z in z]

    plt = meshscatter!(
        ax, pos;
        markersize = sizes,
        marker = Rect3((-0.5, -0.5, eps()), (1, 1, 1)),
        color = z,
        kwargs...
    )

    return plt
end
