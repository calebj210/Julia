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

function onesided_2d(;sizes = [(0,17),(1,5),(4,1)], deg = 6)
    fig = Figure()
    
    for i in 1:length(sizes)
        nodes = onesided_nodes(sizes[i]..., type = BigFloat)
        A = vand(nodes; deg)
        b = zeros(BigFloat, size(A)[1])
        b[3] = 1

        # w = pinv(A) * b
        w = A \ b

        display([w nodes])

        # ax = Axis3(
        #     fig[1,i],
        #     aspect = (1/3,1,1/2),
        #     # yreversed = true,
        #     limits = ((-5.5,5.5), (-0.5,18.5), (0,1.75)),
        #     azimuth = 0.6pi,
        #     elevation = .3,
        #     title = "$(2sizes[i][1] + 1)x$(sizes[i][2] + 1)"
        # )
        #
        # scatter!(
        #     ax, first.(nodes), last.(nodes), abs.(w), 
        #     colorrange = (0, 2), 
        #     color = abs.(w), 
        #     colormap = :viridis,
        # )
        ax = Axis(
            fig[1,i],
            limits = ((-7.5,7.5), (-0.5,30.5)),
            title = "$(2sizes[i][1] + 1)x$(sizes[i][2] + 1)",
            aspect = DataAspect(),
            xticks = -7:7:7,
            xlabel = L"x",
        )

        heatmap!(
            ax, first.(nodes), last.(nodes), abs.(w), 
            colorrange = (1e-3, 200), 
            colormap = :viridis,
            colorscale = log10,
        )
    end

    fig.content[1].ylabel = L"y"

    Colorbar(fig[1,4], limits = (1e-3, 200), colormap = :viridis, scale = log10, label = L"|w_k|")

    rowsize!(fig.layout, 1, Aspect(1,2))
    resize_to_layout!(fig)

    return fig
end
