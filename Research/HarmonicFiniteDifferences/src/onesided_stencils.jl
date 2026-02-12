using CairoMakie, LaTeXStrings
using SparseArrays

include("stencil_solvers.jl")

function onesided_nodes(x::T, y::T, z::T; extended_precision = false, λ = 1) where T <: Integer
    if extended_precision
        return vec([BigFloat.((i,j,λ*k)) for i in -x:x, j in -y:y, k in 0:z])
    else
        return vec([(i,j,λ*k) for i in -x:x, j in -y:y, k in 0:z])
    end
end
onesided_nodes(N; kwargs...) = onesided_nodes(N, N, 2N; kwargs...)

function onesided_sparse_structure(x::Vector, deg)
    A = vand(x, deg)
    B = qr(A')
    idx = abs.(B.R') .> 1e-15
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

onsesided_sparse_structure(N::T, deg) where T<:Integer = onesided_sparse_structure(onesided_nodes(N), deg)
