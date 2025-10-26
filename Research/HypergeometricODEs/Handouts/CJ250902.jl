#=
#   CJ250902 graphics and tests
#
# Author: Caleb Jacobs
# DLM: September 3, 2025
=#

using CairoMakie, ComplexVisuals, LaTeXStrings
using Random, BenchmarkTools
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
    mac = [first(t(a, b, c, z, maclaurin_2f1)) for t ∈ transformations]
    con = [first(t(a, b, c, z, conformal_2f1)) for t ∈ transformations]

    erm, tm = findmin(geterror.(mac, tru))
    erc, tc = findmin(geterror.(con, tru))

    if erm < erc
        er = erm
        mac_or_con = :cross
        tran = tm
    else
        er = erc
        mac_or_con = :circle
        tran = tc
    end

    return (er, mac_or_con, tran)
end

function pert(x)
    if iszero(imag(x))
        return x + rand([-1,1]) * eps(real(x))
    else
        return x + rand([-1,1]) * eps(real(x)) + rand([-1,1]) * eps(imag(x)) * im
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
        marker = vec(getindex.(errsmethod, 2)),
    )

    if isnothing(truerr)
        println("Getting true solution for best of all")
        truerr = johansson_2f1.(a, b, c, zerr; bits = 512)
    end
    errs = besteval.(a, b, c, zerr, truerr)
    surface!(axerror, reim(zerr)..., log10.(first.(errs)))

    return (trumethod, truerr)
end

function test_figures(test, tru_errors = nothing, tru_method = nothing, tru_conds = nothing)
    a, b, c = test

    zmethod = complex_square_grid(2.0, 20)
    zerr = complex_square_grid(2.0, 50)

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

    return (fig, (tru_error, tru_method, tru_conds))
end

function histograms(;N = 10000, arng = 25, brng = 25, crng = 25, zrng = 3, seed = 997)
    Random.seed!(seed)

    # Setup random tests
    as = arng * (1 .- 2rand(N))
    bs = brng * (1 .- 2rand(N))
    cs = crng * (1 .- 2rand(N))
    tests = Vector{Tuple{Float64, Float64, Float64, ComplexF64}}()
    zs = zrng * (2(rand(N) + rand(N) * im) .- (1 + im))

    print("Getting tests: ")
    tru = Vector{ComplexF64}()
    alt = Vector{Float64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = johansson_2f1(a,b,c,z; bits = 512)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))

        altval = johansson_2f1(pert.((a, b, c, z))...; bits = 512)
        push!(alt, geterror(altval, val))
        
    end
    println("count: $(length(tests))")

    # Generate graphics and results
    set_theme!(theme_latexfonts())
    fig = Figure()
    print("Running best of all: ")
    err = [besteval(test..., tr)[1] for (test, tr) ∈ zip(tests, tru)]
    # joerr = geterror.([johansson_2f1(test...; bits = 53) for test ∈ tests], tru)
    println("done!")
    print("Making graphics: ")

    ax1 = Axis(fig[1,1], 
        limits = (nothing, (0,1)), 
        xscale = log10, 
        xlabel = "Relative Error", 
        title = "Best Series/Conformal Error",
    )
    ax2 = Axis(fig[1,2], 
        limits = (nothing, (0,1)), 
        xscale = log10, 
        xlabel = "Relative Difference", 
        # xlabel = "Relative Error", 
        title = "Perturbed 2F1",
    )

    bin = 10.0 .^ (-16:2:2)

    hist!(ax1, err, 
        bins = bin, 
        # color = :values,
        color = :black,
        normalization = :probability
    )
    # hist!(ax2, joerr, 
    hist!(ax2, alt, 
        bins = bin, 
        # color = :values, 
        color = :black,
        normalization = :probability
    )

    rowsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)
    println("done!")

    return fig
end
