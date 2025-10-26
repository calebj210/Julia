#=
#   Comparison tests for the conformal mapping and Maclaurin series.
#
# Author: Caleb Jacobs
# DLM: October 22, 2025
=#

include("../Base2F1.jl")
include("../pFq.jl")

using GLMakie, LaTeXStrings
using ComplexVisuals

function geterr(f, t)
    val = abs(f - t) / (abs(t) + eps(abs(t)))

    if isinf(val) || isnan(val) || val > 1
        return 1e0
    elseif val < eps()
        return eps()
    else
        return val
    end
end

function compare(a, b, c, z)
    con = conformal_2f1.(a, b, c, z)
    mac = maclaurin_2f1.(a, b, c, z)
    tru = johansson_2f1.(a, b, c, z; bits = 512)

    erc = geterr.(con, tru)
    erm = geterr.(mac, tru)
    best = mapslices(x -> findmin(x)[2], cat(erc, erm; dims = 3), dims = 3)[:,:,1]

    fig = Figure()
    axc = Axis3(
        fig[1,1],
        title = "Conformal Mapping",
        xlabel = "Re(z)",
        ylabel = "Im(z)",
    )
    axm = Axis3(
        fig[1,2],
        title = "Maclaurin Series",
        xlabel = "Re(z)",
        ylabel = "Im(z)",
    )
    axb = Axis(
        fig[2,1:2],
        title = "Best val",
        xlabel = "Re(z)",
        ylabel = "Im(z)",
    )

    surface!(axc, reim(z)..., log10.(erc), colorrange = (log10(eps()), 1))
    surface!(axm, reim(z)..., log10.(erm), colorrange = (log10(eps()), 1))
    plt = heatmap!(axb, z.real, z.imag, best)

    Colorbar(fig[2,3], plt)

    resize_to_layout!(fig)

    return fig
end
