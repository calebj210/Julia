#=
# Tests for inspected lower accuracy random tests
#
# Author: Caleb Jacobs
# DLM: July 21, 2025
=#

using ComplexVisuals, CairoMakie, LaTeXStrings
using Random
include("pFq.jl")

plog(z) = sign(z) * log10(abs(z) + 1)

function grid_errors(test)
    a, b, c = test[1:3]
    z = complex_square_grid(3,300)

    println("Getting true solution")
    tru = johansson_2f1.(a, b, c, z, bits = 512)
    println("Getting Taylor solution")
    val = taylor_2f1.(a, b, c, z)

    err = clean_error.(val, tru)

    set_theme!(complex_theme)
    fig = Figure()
    ax1 = Axis3(fig[1,1],
        title = "Relative Error",
        zlabel = L"\log_{10}|\text{Error}|",
        zlabeloffset = 40,
        viewmode = :fit
    )
    ax2 = Axis3(fig[1,2],
        title = "Phase Portrait",
        zlabel = L"pseudolog$_{10}(f)$",
        zlabeloffset = 40,
    )
    surface!(ax1, reim(z)..., log10.(err))
    complexsurface!(ax2, z, sign.(tru) .* log10.(abs.(tru) .+ 1))

    axiscolorwheel(ax2, position = :rt)

    rowsize!(fig.layout, 1, Aspect(1, 1))

    resize_to_layout!(fig)

    return fig
end

function curious_tests()
    tests = (
        (-7.004297927757014, 12.140335057326817, -19.501876376538124),  # z = 1 root wall   (out)
        (17.708272572407935, -7.277152107620976, -24.056815489580536),  # z = 1 root wall   (out)
        (4.814400989782125, -16.4025275971892, -22.299930456581336),    # Root wall         (out)
        (-6.286902255001703, 1.7926091191509619, -14.753432765993685),  # Mild root walls   (out)
        (0.3413776181071859, 24.169612904438516, -13.54665644705042),   # Root wall         (in)
        (-2.1909384761924224, -16.27009747370081, -8.918934458011163),  # Root wall         (in)
        (0.6545186879146792, 22.010065823632917, -20.67432189526394),   # Root wall         (in)
        (-2.320440016143699, -23.526716481523753, 16.37180276700314),   # Direction and wall(in)
        (-21.50537648336128, -1.759360586651454, -15.922872131260618),  # Root wall         (in)
        (17.708272572407935, -7.277152107620976, -24.056815489580536),  # z = 1 init        (in)
    )

    return tests
end

function random_test(;N = 10000, arng = 25, brng = 25, crng = 25, zrng = 2, seed = 997, complextest = false)
    Random.seed!(seed)

    # Setup random tests
    if complextest
        as = arng * complexrand(N)
        bs = brng * complexrand(N)
        cs = crng * complexrand(N)
        tests = Vector{NTuple{4, ComplexF64}}()
    else
        as = arng * (1 .- 2rand(N))
        bs = brng * (1 .- 2rand(N))
        cs = crng * (1 .- 2rand(N))
        tests = Vector{Tuple{Float64, Float64, Float64, ComplexF64}}()
    end
    zs = zrng * complexrand(N)

    print("Getting tests: ")
    tru = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = johansson_2f1(a,b,c,z; bits = 512)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))
    end
    println("count: $(length(tests))")

    # Generate graphics and results
    set_theme!(theme_latexfonts())
    fig = Figure()
    err = clean_error.([taylor_2f1(test...) for test ∈ tests], tru)

    ax = Axis(fig[1,1], 
        limits = (nothing, (0,1)), 
        xscale = log10, 
        xlabel = "Relative Error", 
    )

    # bin = 10.0 .^ (-16:2:2)
    bin = 10.0 .^ (-16:1:0)

    hist!(ax, err, 
        bins = bin, 
        color = :values, 
        normalization = :probability
    )

    sorted_tests = sort_errors(tests, err, 10.0 .^ (-15:1:0))

    return fig, sorted_tests
end

function sort_errors(tests, errs, bins)
    sorted_tests = ntuple(x -> Vector{Tuple{Float64, Float64, Float64, ComplexF64}}(), length(bins))

    for (err, test) ∈ zip(errs, tests)
        for (i, bin) ∈ enumerate(bins)
            if err <= bin
                push!(sorted_tests[i], test)
                break
            end
        end
    end

    # return sorted_tests
    return Dict(zip(round.(Int, log10.(bins)), sorted_tests))
end

function clean_error(f,t)
    err = abs(f - t) / (abs(t) + eps())

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-16
        err = 1e-16
    end

    return err
end

function complexrand(N)
    vals = -(1 + im) .+ 2rand(ComplexF64, N)

    return vals
end

function split_transformation(a, b, c)
    z = complex_square_grid(3, 300)

    gl = gamma(c - a - b) * gamma(c) / gamma(c - a) / gamma(c - b)
    gr = gamma(a + b - c) * gamma(c) / gamma(a) / gamma(b) * (1 .- z).^(c - a - b)

    fl = taylor_2f1.(a, b,         a + b - c + 1, 1 .- z)
    fm = taylor_2f1.(a, b, c, z)
    fr = taylor_2f1.(c - a, c - b, c - a - b + 1, 1 .- z)
    

    fig = Figure()
    axl = Axis3(fig[1,1], title = "Left Full")
    axm = Axis3(fig[1,2], title = "Combined Full")
    axr = Axis3(fig[1,3], title = "Right Full")

    axbl = Axis3(fig[2,1], title = "Left")
    axbm = Axis3(fig[2,2], title = "Computed Full")
    axbr = Axis3(fig[2,3], title = "Right")

    complexsurface!(axl, z, plog.(gl .* fl))
    complexsurface!(axm, z, plog.(gl .* fl + gr .* fr))
    complexsurface!(axr, z, plog.(gr .* fr))

    complexsurface!(axbl, z, plog.(fl))
    complexsurface!(axbm, z, plog.(fm))
    complexsurface!(axbr, z, plog.(fr))

    return fig
end

function conditioning(a, b, c)
    z = complex_square_grid(3,300)

    α, β, γ = (a, b, c) .+ eps.((a, b, c))

    f = johansson_2f1.(a, b, c, z, bits = 512)
    ϕ = johansson_2f1.(α, β, γ, z, bits = 512)

    dif = clean_error.(ϕ, f)

    fig = Figure()
    axl = Axis3(fig[1,1],
        title = "Phase Portrait",
        zlabel = L"pseudolog$_{10}$(F)",
        zlabeloffset = 40,
        viewmode = :fit
    )
    axr = Axis3(fig[1,2],
        title = "Relative Difference",
                zlabel = L"\log_{10}|F - F̃| / |F|",
        zlabeloffset = 40,
        viewmode = :fit
    )

    complexsurface!(axl, z, plog.(f))
    surface!(axr, reim(z)..., log10.(dif))

    axiscolorwheel(axl, position = :rt)

    rowsize!(fig.layout, 1, Aspect(1, 1))
    resize_to_layout!(fig)

    return fig
end
