#=
#   Functions for generating graphics in the paper
#
# Author: Caleb Jacobs
# DLM: November 6, 2025
=#

using CairoMakie, ComplexVisuals, LaTeXStrings
using Random, BenchmarkTools

include("../Comparison2F1.jl")
include("../pFq.jl")

# 2F1 methods and names
funcs = (comparison_2f1, weniger_2f1, (a,b,c,z) -> johansson_2f1(a, b, c, z, bits = 53), mathematica_2f1, matlab_2f1)
names = ("Conformal", "Levin-Type", "Johansson", "Mathematica", "MATLAB")
indices = ((3,1), (1,1), (1,2), (2,1), (2,2))

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

    for _ ∈ 1 : N
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
    ecrng = (-16,-8)
    etickrng = -15:3:-9
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
        zlabel = L"\text{abs}(f)",
        title = L"\text{Abs-Arg of } {_2}F_1(1,-9/2;-9/4;z)",
        elevation = 0.2π,
        xticklabelsize = 10, yticklabelsize = 10, zticklabelsize = 10,
        xlabeloffset = 20, ylabeloffset = 20, zlabeloffset = 53,

    )
    complexsurface!(ax, z, tru)
    axiscolorwheel(ax, width = Relative(.15), height = Relative(.15), halign = .85)

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
    alt = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = johansson_2f1(a,b,c,z; bits = 512)
        valalt = johansson_2f1(perturb.((a,b,c,z))...; bits = 1024)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))
        push!(alt, convert(ComplexF64, valalt))
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
            color = :black, 
            normalization = :probability
        )
    end
    ax = Axis(fig[3,2], 
        limits = (nothing, (0,1)), 
        xscale = log10, 
        xlabel = "Relative Difference", 
        title = "Sensitivity",
    )

    bin = 10.0 .^ (-16:2:2)

    dif = clean_error.(alt, tru)
    hist!(ax, dif, 
        bins = bin, 
        color = :black, 
        normalization = :probability
    )

    return fig
end

# Alternate branch
function alternate_branch(a = 1.1, b = .5, c = 1.2)
    z = ComplexGrid(range(-1,2,301), range(-1.5,1.5,301))

    f1 = comparison_2f1.(a, b, c, z)
    f2 = [π * gamma(c) / sinpi(b - a) * (
        cispi((imag(z) < 0 ? 1 : -1) * a) * z^(-a) / gamma(b) / gamma(c - a) / gamma(a - b + 1) *
        comparison_2f1(a, a - c + 1, a - b + 1, 1 / z) -
        cispi((imag(z) < 0 ? 1 : -1) * b) * z^(-b) / gamma(a) / gamma(c - b) / gamma(b - a + 1) *
        comparison_2f1(b, b - c + 1, b - a + 1, 1 / z)
    ) for z ∈ z]

    # f1[isnan.(f1)] .= 0
    # f2[isnan.(f2)] .= 0

    relims = (-2,2)
    imlims = (-2,2)

    # f1[(real(z) .>= 1) .&& (imag(z) .== 0)] .= NaN + NaN * im
    # f2[(imag(z) .== 0)] .= NaN + NaN * im

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax1 = Axis3(fig[1,1],
        title = "Real Part",
        titlesize = 25,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"\text{abs}($f$)",
        xticklabelsize = 15, yticklabelsize = 15, zticklabelsize = 15,
        xlabelsize = 20, ylabelsize = 20, zlabelsize = 20,
        xlabeloffset = 20, ylabeloffset = 20, zlabeloffset = 25,
        limits = (nothing, nothing, relims),
        azimuth = .275π,
        elevation = .15π,
    )
    ax2 = Axis3(fig[1,2],
        title = "Imaginary Part",
        titlesize = 25,
        xlabel = L"\mathrm{Re}(z)", ylabel = L"\mathrm{Im}(z)", zlabel = L"\text{abs}($f$)",
        xticklabelsize = 15, yticklabelsize = 15, zticklabelsize = 15,
        xlabelsize = 20, ylabelsize = 20, zlabelsize = 20,
        xlabeloffset = 20, ylabeloffset = 20, zlabeloffset = 25,
        limits = (nothing, nothing, imlims),
        azimuth = .275π,
        elevation = .15π,
    )

    surface!(ax1, reim(z[:,1:150])..., real(f1[:,1:150]), alpha = .9, colormap = :blues, colorrange = relims, label = "Principle Sheet")
    surface!(ax1, reim(z[:,151:end])..., real(f1[:,151:end]), alpha = .9, colormap = :blues, colorrange = relims)
    surface!(ax1, reim(z[:,1:150])..., real(f2[:,1:150]), alpha = .9, colormap = :reds, colorrange = relims, label = "Second Sheet")
    surface!(ax1, reim(z[:,151:end])..., real(f2[:,151:end]), alpha = .9, colormap = :reds, colorrange = relims)
    surface!(ax2, reim(z[:,1:150])..., imag(f1[:,1:150]), alpha = .9, colormap = :blues, colorrange = imlims)
    surface!(ax2, reim(z[:,151:end])..., imag(f1[:,151:end]), alpha = .9, colormap = :blues, colorrange = imlims)
    surface!(ax2, reim(z[:,1:150])..., imag(f2[:,1:150]), alpha = .9, colormap = :reds, colorrange = imlims)
    surface!(ax2, reim(z[:,151:end])..., imag(f2[:,151:end]), alpha = .9, colormap = :reds, colorrange = imlims)

    lines!(ax1, reim(z[:,150])..., real(f1[:,150]), color = :black, linewidth = 3)
    lines!(ax1, reim(z[:,150])..., real(f2[:,150]), color = :black, linewidth = 3)
    lines!(ax1, reim(z[:,151])..., real(f2[:,151]), color = :black, linewidth = 3)
    #
    lines!(ax2, reim(z[:,150])..., imag(f1[:,150]), color = :black, linewidth = 3)
    lines!(ax2, reim(z[:,150])..., imag(f2[:,150]), color = :black, linewidth = 3)
    lines!(ax2, reim(z[:,151])..., imag(f2[:,151]), color = :black, linewidth = 3)

    axislegend(ax2, ax1, position = :rt, labelsize = 20)

    rowsize!(fig.layout, 1, Aspect(1,1))

    resize_to_layout!(fig)

    return fig
end

# Sensitivity plot
perturb(a::Float64) = a + eps(a)
perturb(a::ComplexF64) = perturb(real(a)) + im * perturb(imag(a))
function sensitivity(;N = 10000, arng = 25, brng = 25, crng = 25, zrng = 2, seed = 997, complextest = false)
    Random.seed!(seed)

    # Setup random tests
    if complextest
        as = arng * complexrand(N)
        bs = brng * complexrand(N)
        cs = crng * complexrand(N)
    else
        as = arng * (1 .- 2rand(N))
        bs = brng * (1 .- 2rand(N))
        cs = crng * (1 .- 2rand(N))
    end
    zs = zrng * complexrand(N)

    print("Getting tests: ")
    tru = Vector{ComplexF64}()
    alt = Vector{ComplexF64}()
    comps = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = johansson_2f1(a,b,c,z; bits = 1024)
        valalt = johansson_2f1(perturb.((a,b,c,z))...; bits = 1024)
        comp = comparison_2f1(a,b,c,z)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tru, convert(ComplexF64, val))
        push!(alt, convert(ComplexF64, valalt))
        push!(comps, comp)
    end

    # Generate graphics and results
    set_theme!(theme_latexfonts())
    fig = Figure()
    err = clean_error.(alt, tru)
    myerr = clean_error.(comps, tru)
    bin = 10.0 .^ (-16:2:2)
    ax = Axis(fig[1,1], 
              limits = (nothing, (0,1)), 
              xscale = log10, 
              xlabel = "Relative Difference", 
              title = "Sensitivity",
    )
    hist!(ax, err, 
          bins = bin, 
          color = :black, 
          normalization = :probability
    )
    ax = Axis(fig[1,2], 
              limits = (nothing, (0,1)), 
              xscale = log10, 
              xlabel = "Relative Error", 
              title = "Conformal",
    )
    hist!(ax, myerr, 
          bins = bin, 
          color = :black, 
          normalization = :probability
    )

    rowsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end
