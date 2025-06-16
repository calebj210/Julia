#=
# Makie visuals for 2D harmonic basis functions
# 
# Author: Caleb Jacobs
# DLM: June 15, 2025
=#

using GLMakie, LaTeXStrings
include("2D_Basis.jl")

function plotbasisfunction(N)
    x = range(-1, 1, 100)
    y = range(-1, 1, 100)

    fig = Figure()
    ax = Axis3(fig[1,1],
        xlabel = "x",
        ylabel = "y",
        zlabel = "f(x,y)",
    )
    surface!(ax, x, y, (x,y) -> harmonic(N, x, y))

    return fig
end

function plotcardinalfunction(N, n = N^2 ÷ 2 + 1)
    x = range(-N ÷ 2, N ÷ 2, 300)
    w = cardinalweights(N)[:, n]

    fig = Figure()
    ax = Axis3(fig[1,1],
        xlabel = "x",
        ylabel = "y",
        zlabel = "f(x,y)",
    )
    surface!(ax, x, x, (x, y) -> log10.(abs(cardinal(N, w, x, y))))

    return fig
end

function plotcardinalfunctions(N; shownodes = false, showwireframe = false, scale = identity)
    if showwireframe
        x = range(-N ÷ 2, N ÷ 2, 25)
    else
        x = range(-N ÷ 2, N ÷ 2, 100)
    end

    w = cardinalweights(N)
    if shownodes
        nodes = gridnodes(N)
    end

    fig = Figure()
    for i ∈ 1:N, j ∈ 1:N
        card(x, y) = scale(cardinal(N, w[:, i + (j - 1) * N], x, y))

        ax = Axis3(fig[i,j],
            xlabel = "x",
            ylabel = "y",
            zlabel = "f(x,y)",
        )

        if showwireframe
            X = collect(x)
            Z = [card.(x, y) for x ∈ x, y ∈ x]

            wireframe!(ax, X, X, Z)
        else
            surface!(ax, x, x, card)
        end

        if shownodes
            colors = repeat([:green], N^2)
            colors[i + (j - 1) * N] = :red
            scatter!(ax, 
                first.(nodes), last.(nodes), 
                card.(first.(nodes), last.(nodes)), 
                color = colors,
                markersize = 15,
            )
        end
    end

    for i ∈ 1:N
        colsize!(fig.layout, i, Aspect(1, 1))
    end
    resize_to_layout!(fig)

    return fig
end

function plotcardinalfunctionzeros(N; colorcontours = false, shownodes = true)
    x = range(-N ÷ 2, N ÷ 2, 100)
    w = cardinalweights(N)
    if shownodes
        nodes = gridnodes(N)
    end

    fig = Figure()
    for i ∈ 1:N, j ∈ 1:N
        ax = Axis(fig[i,j])
        hidedecorations!(ax)
        if colorcontours
            contourf!(ax, x, x, (x,y) -> Makie.pseudolog10(cardinal(N, w[:, i + (j - 1) * N], x, y)),
                levels = -20:2:20,
            )
        else
            contour!(ax, x, x, (x,y) -> cardinal(N, w[:, i + (j - 1) * N], x, y),
                levels = [0],
                color = :black,
            )
        end
        if shownodes
            colors = repeat([:green], N^2)
            colors[i + (j - 1) * N] = :red
            scatter!(ax, first.(nodes), last.(nodes), color = colors)
        end
    end
    if colorcontours
        Colorbar(fig[N ÷ 2:N ÷ 2 + 2, N + 1], 
            limits = (-20,20), 
            size = 35,
            colormap = cgrad(:viridis, 20, categorical = true),
            ticks = (-20:10:20, [p == 0 ? "0" : p < 0 ? 
                latexstring("1 - 10^{", -p, "}") : latexstring("10^{", p, "} - 1") 
                for p ∈ -20:10:20]
            ),
            ticklabelsize = 25,
        )
    end

    for i ∈ 1:N
        colsize!(fig.layout, i, Aspect(1, 1))
    end
    resize_to_layout!(fig)

    return fig
end
