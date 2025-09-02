#=
#   CJ250902 graphics and tests
#
# Author: Caleb Jacobs
# DLM: September 2, 2025
=#

using CairoMakie, ComplexVisuals, LaTeXStrings
include("../pFq.jl")
include("../Conformal2F1.jl")

legend_colors = cgrad(:Set1_8, 8, categorical = true)

function geterror(f, g)
    err = abs(f - g) / (abs(g) + eps(abs(f)))
    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err < 1e-16
        err = 1e-16
    end
    return err
end

function besteval(a, b, c, z, tru)
    mac = [first(t.(a, b, c, z, maclaurin_2f1)) for t ∈ transformations]
    con = [first(t.(a, b, c, z, conformal_2f1)) for t ∈ transformations]

    erm, tm = findmin(geterror.(mac, tru))
    erc, tc = findmin(geterror.(con, tru))

    if erm < erc
        er = erm
        tran = tm
        mac_or_con = 1
    else
        er = erc
        tran = tc
        mac_or_con = 2
    end

    return (er, mac_or_con, tran)
end

function pert(x)
    if iszero(imag(x))
        return x + (2rand() - 1) * eps(real(x))
    else
        return x + (2rand() - 1) * eps(real(x)) + (2rand() - 1) * eps(imag(x)) * im
    end
end

function condition_plot!(ax, a, b, c, z; tru = nothing)
    if isnothing(tru)
        println("Generating condition plot...")
        f = johansson_2f1.(a, b, c, z; bits = 512)
        g = johansson_2f1.(pert(a), pert(b), pert(c), pert.(z); bits = 512)

        tru = geterror.(g, f)
    end

    surface!(ax, reim(z)..., log10.(tru))

    return tru
end

function method_plot!(axerror, axmethod, a, b, c, zmethod, zerr; trumethod = nothing, truerr = nothing)
    if isnothing(trumethod)
        println("Getting true solution for best of all")
        trumethod = johansson_2f1.(a, b, c, zmethod; bits = 512)
    end

    println("Getting best evaluations")
    errsmethod = besteval.(a, b, c, zmethod, trumethod)

    scatter!(axmethod, reim(vec(zmethod))...,
        color = vec(last.(errsmethod)),
        colormap = :Set1_8,
        colorrange = (1,8),
        marker = vec([idx == 1 ? :cross : :circle for idx ∈ getindex.(errsmethod, 2)]),
    )

    if isnothing(truerr)
        println("Getting true solution for best of all")
        truerr = johansson_2f1.(a, b, c, zerr; bits = 512)
    end
    errs = besteval.(a, b, c, zerr, truerr)
    surface!(axerror, reim(zerr)..., log10.(first.(errs)))

    return (trumethod, truerr)
end

function histograms()

    return nothing
end

function test_figures(test, tru_errors = nothing, tru_method = nothing, tru_conds = nothing)
    a, b, c = test

    zmethod = complex_square_grid(3.0, 20)
    zerr = complex_square_grid(3.0, 300)

    set_theme!(theme_latexfonts())
    fig = Figure()
    errors = Axis3(fig[1,1],
        title = "Relative Error",
        zlabel = L"\log_{10}|\text{Error}|",
        zlabeloffset = 40,
        viewmode = :fit
    )
    conds = Axis3(fig[1,2],
        title = "Conditioning",
        zlabel = L"\log_{10}|\text{Error}|",
        zlabeloffset = 40,
        viewmode = :fit
    )
    methods = Axis(fig[2,1],
        aspect = DataAspect(),
        title = "Method Used",
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
    )

    tru_method, tru_error = method_plot!(errors, methods, a, b, c, zmethod, zerr; trumethod = tru_method, truerr = tru_errors)
    tru_conds = condition_plot!(conds, a, b, c, zerr; tru = tru_conds)

    Legend(
        fig[2,2],
        [
            MarkerElement(color = :black, marker = :cross),
            MarkerElement(color = :black, marker = :circle),
            [MarkerElement(marker = :rect, color = legend_colors[i]) for i ∈ 1:8]...
        ],
        [
            "Maclaurin Series",
            "Conformal Mapping",

            L"z",
            L"z \mapsto z",
            L"$z \mapsto z / (z - 1)$ (a)",
            L"$z \mapsto z / (z - 1)$ (b)",
            L"z \mapsto 1 - z",
            L"z \mapsto 1 - 1 / z",
            L"z \mapsto 1 / z",
            L"z \mapsto 1 / (1 - z)"
        ],
        rowgap = -3,
    )

    resize_to_layout!(fig)

    return (fig, tru_error, tru_method, tru_conds)
end
