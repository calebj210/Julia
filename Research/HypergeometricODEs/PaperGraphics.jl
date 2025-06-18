#=
#   Functions for generating graphics in the paper
#
# Author: Caleb Jacobs
# DLM: June 18, 2025
=#

using CairoMakie, ComplexVisuals, LaTeXStrings
using Random, BenchmarkTools
include("pFq.jl")
include("PaperTests/Slevinsky2F1.jl")

# 2F1 methods and names
const funcs = (taylor_2f1, weniger_2f1, (a,b,c,z) -> johansson_2f1(a, b, c, z, bits = 53), mathematica_2f1)
const names = ("Taylor", "Levin-Type", "Johansson", "Mathematica")
const indices = ((1,1), (1,2), (2,1), (2,2))

# Helper functions
function clean_error(f,t)
    err = abs.((f - t) ./ t)

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-16
        err = 1e-16
    end

    return err
end

function average_time(a, b, c, z, f, N = 5; microseconds = true)
    time = 0.0

    for n ∈ 1 : N
        time += @elapsed f(a, b, c, z)
    end

    time /= N

    if microseconds
        time = round(Int, time * 1e6) # Convert to microseconds
    end

    return time
end

function complexrand(N)
    vals = -(1 + im) .+ 2rand(ComplexF64, N)

    return vals
end

## Tests
# Slevinsky grid test: 2F1(1,-9/2;-9/4;z)
function slevinsky_grid_test(tru = nothing)
    z = complex_square_grid(5, 300)
    test = (1, -9/2, -9/4, z)

    if isnothing(tru)
        println("Getting true solution")
        tru = johansson_2f1.(test...)
    end

    # Generate graphics
    set_theme!(theme_latexfonts())
    ecrng = (-16,-5)
    etickrng = -15:3:-6
    eticks = (etickrng, [latexstring("10^{", p, "}") for p ∈ etickrng])

    tcrng = (-6,-2)
    ttickrng = -5:-3
    tticks = (ttickrng, [latexstring("10^{", p, "}") for p ∈ ttickrng])

    default_settings = (;
        zticks = eticks,
        zlabelsize = 12,
        limits = (nothing, nothing, ecrng),
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        xticklabelsize = 10, yticklabelsize = 10, zticklabelsize = 10,
        xlabeloffset = 20, ylabeloffset = 20, zlabeloffset = 40,
    )

    errors  = Figure()
    timings = Figure()

    # Surfaces
    for (f, name, idx) ∈ zip(funcs, names, indices)
        println("Running ", name)
        println("\t Errors")
        err = clean_error.(f.(test...), tru)
        println("\t Timings")
        time = average_time.(test..., f, microseconds = false)

        axe = Axis3(errors[idx...];
            default_settings...,
            title = string(name, " Errors"),
        )
        axt = Axis3(timings[idx...];
            default_settings...,
            limits = (nothing, nothing, tcrng),
            zticks = tticks,
            zlabel = "Time",
            title = string(name, " Timings"),
        )

        surface!(axe, reim(z)..., log10.(err), colorrange = ecrng)
        surface!(axt, reim(z)..., log10.(time), colorrange = tcrng)
    end

    # Colorbar
    Colorbar(errors[:,3], limits = ecrng, ticks = eticks)
    Colorbar(timings[:,3], limits = tcrng, ticks = tticks)

    # Formatting
    resize_to_layout!(errors)
    resize_to_layout!(timings)

    # Phase portrait
    ph = Figure()
    ax = Axis3(ph[1,1];
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"|f|",
        elevation = 0.2π,
        xticklabelsize = 10, yticklabelsize = 10, zticklabelsize = 10,
        xlabeloffset = 20, ylabeloffset = 20, zlabeloffset = 40,

    )
    complexsurface!(ax, z, tru)

    return (;errors, timings, ph, tru)
end

# Random histograms
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
        val = arb_2f1(ArbComplex.((a,b,c,z), bits = 512)...)
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
    for (f, name, idx) ∈ zip(funcs, names, indices)
        println("Running ", name)
        println("\tErrors")
        err = clean_error.([f(test...) for test ∈ tests], tru)
        print("\tTimings: ")
        time = [average_time(test..., f) for test ∈ tests]
            println("mean(", round(mean(time)), "), median(", round(median(time)), ")")

        ax = Axis(fig[idx...], 
            limits = (nothing, (0,1)), 
            xscale = log10, 
            xlabel = "Relative Error", 
            title = name,
        )

        bin = 10.0 .^ (-16:2:2)

        hist!(ax, err, 
            bins = bin, 
            color = :values, 
            normalization = :probability
        )
    end

    return fig
end
