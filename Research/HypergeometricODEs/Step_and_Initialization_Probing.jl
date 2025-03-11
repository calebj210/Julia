#= 
# Tests for inspecting failures of the Taylor method
#
# Author: Caleb Jacobs
# DLM: March 11, 2025
=#

include("pFq.jl")
include("PaperTests/Slevinsky2F1.jl")

using CairoMakie, ComplexVisuals, LaTeXStrings
import Random.seed!

function init_step_error(a, b, c, z)
    init_h = 10 .^ range(-2, log10(1), 300)
    step_h = 10 .^ range(-3, log10(10), 300)

    vals = [taylor_2f1(a, b, c, z; init_max = ih, step_max = sh) for ih ∈ init_h, sh ∈ step_h]
    tru = johansson_2f1(a, b, c, z, bits = 2024)

    errs = abs.((vals .- tru) ./ tru)
    errs[errs .> 1] .= 1
    errs[errs .<= 1e-17] .= 1e-17

    fig = Figure()
    ax  = Axis(fig[1,1], 
               xlabel = "Initial Step", 
               ylabel = "Max Step",
               title  = "Initialization and Step Size Error",
               xscale = log10,
               yscale = log10)
    plt = heatmap!(ax, init_h, step_h, errs, 
                   colorscale = log10,
                   )
    lines!(ax, init_h, init_h, color = :orange)

    Colorbar(fig[1,2], plt)
    colsize!(fig.layout, 1, Aspect(1, 1))
    resize_to_layout!(fig)

    return fig
end

function global_error(a, b, c; z = nothing, tru = nothing)
    if isnothing(z)
        z = complex_square_grid(5, 300)
    end

    if isnothing(tru)
        print("Computing true solution...")
        tru = johansson_2f1.(a, b, c, z, bits = 1024)
        println("done!")
    end

    print("Computing Taylor...")
        # val1 = taylor_2f1.(a, b, c, z, two_step = false)
        # val2 = taylor_2f1.(a, b, c, z, two_step = true)
        val1 = _2f1.(a, b, c, z, two_step = false)
        val2 = _2f1.(a, b, c, z, two_step = true)
    println("done!")

    errs1 = clean_error.(val1, tru)
    errs2 = clean_error.(val2, tru)

    set_theme!(theme_latexfonts())
    fig = Figure()

    errs = GridLayout(fig[1,1])
    phases = GridLayout(fig[2,1])
    ax11 = Axis(errs[1,1], xautolimitmargin = (0,0), yautolimitmargin = (0,0), title = "Single Step Error", ylabel = L"Im(z)")
    ax12 = Axis(errs[1,2], xautolimitmargin = (0,0), yautolimitmargin = (0,0), title = "Two Step Error")
    ax21 = Axis(phases[1,1], xautolimitmargin = (0,0), yautolimitmargin = (0,0), title = "Single Step Phase", xlabel = L"Re(z)", ylabel = L"Im(z)")
    ax22 = Axis(phases[1,2], xautolimitmargin = (0,0), yautolimitmargin = (0,0), title = "Two Step Phase", xlabel = L"Re(z)")
    ax3 = Axis(phases[1,3], xautolimitmargin = (0,0), yautolimitmargin = (0,0), title = "True Phase", xlabel = L"Re(z)")

    plt11 = heatmap!(ax11, z, errs1, colorscale = log10, colorrange = (1e-17,1))
    plt12 = heatmap!(ax12, z, errs2, colorscale = log10, colorrange = (1e-17,1))
    plt21 = phase!(ax21, z, val1)
    plt22 = phase!(ax22, z, val2)
    plt3 = phase!(ax3, z, tru)

    complex_color_wheel!(ax21.scene)
    complex_color_wheel!(ax22.scene)
    complex_color_wheel!(ax3.scene)

    colsize!(phases, 1, Aspect(1,1))
    colsize!(phases, 2, Aspect(1,1))
    colsize!(phases, 3, Aspect(1,1))
    colsize!(errs, 1, Aspect(1,1))
    colsize!(errs, 2, Aspect(1,1))

    Colorbar(errs[1,3], plt11)

    resize_to_layout!(fig)

    return (fig, z, tru)
end

function random_failed_tests(a = 0, b = 0, c = 0, z = 0; N = 10000, arng = 30, brng = 30, crng = 30, zrng = 1, seed = 997, two_step = true)
    seed!(seed)

    # Setup random tests
    # as = a .+ arng * (1 .- rand(N))
    # bs = b .+ brng * (1 .- rand(N))
    # cs = c .+ crng * (1 .- rand(N))
    # zs = z .+ zrng * (1 .- rand(N))
    as = a .+ arng * complexrand(N)
    bs = b .+ brng * complexrand(N)
    cs = c .+ crng * complexrand(N)
    zs = z .+ zrng * complexrand(N)

    print("Getting tests... ")
    tests = Vector{NTuple{4, ComplexF64}}()
    tru = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = arb_2f1(ArbComplex.((a,b,c,z), bits = 512)...)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))
    end
    println("done\nTest Count: $(length(tests))")

    # Evaluate each test for accuracy 
    println("\nRunning accuracy tests:")
    print("\tTaylor: ")
    ta = [taylor_2f1(test..., two_step = two_step)               for test ∈ tests]
    tae = clean_error.(ta, tru)
    taf = tests[isone.(tae)]
    print("done\n\tTransforms: ")
    tr = [_2f1(test..., two_step = two_step)                     for test ∈ tests]
    tre = clean_error.(tr, tru)
    print("done\n\tLevin: ")
    le = [weniger_2f1(test...)              for test ∈ tests]
    lee = clean_error.(le, tru)
    print("done\n\tJohansson: ")
    jo = [johansson_2f1(test..., bits = 53) for test ∈ tests]
    joe = clean_error.(jo, tru)
    # print("done\n\tMathematica: ")
    # ma = [mathematica_2f1(test...)          for test ∈ tests]
    # mae = clean_error.(ma, tru)
    # maf = tests[isone.(mae)]
    println("done")
    
    print("Collecting tests: ")
    fls = taf
    println("done")

    # Generate histograms
    fig = Figure()
    top = GridLayout(fig[1,1])
    bot = GridLayout(fig[2,1])
    axt = Axis(top[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Taylor")
    axr = Axis(top[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Taylor Transformations")
    axl = Axis(bot[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Levin-Type")
    axj = Axis(bot[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Johansson")
    # axm = Axis(bot[1,3], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Mathematica")

    bin = 10.0 .^ (-17:1)

    hist!(axt, tae, bins = bin, color = :values, normalization = :probability)
    hist!(axr, tre, bins = bin, color = :values, normalization = :probability)
    hist!(axl, lee, bins = bin, color = :values, normalization = :probability)
    hist!(axj, joe, bins = bin, color = :values, normalization = :probability)
    # hist!(axm, mae, bins = bin, color = :values, normalization = :probability)
    
    rowsize!(fig.layout, 1, Auto(1))
    rowsize!(fig.layout, 2, Auto(1))

    # return (; fig, fails = fails.taf)
    return (; fig, fls)
end

function random_histograms(a = 0, b = 0, c = 0, z = 0; N = 10000, arng = 30, brng = 30, crng = 30, zrng = 1, seed = 997)
    seed!(seed)

    # Setup random tests
    # as = a .+ arng * (1 .- rand(N))
    # bs = b .+ brng * (1 .- rand(N))
    # cs = c .+ crng * (1 .- rand(N))
    # zs = z .+ zrng * (1 .- rand(N))
    as = a .+ arng * complexrand(N)
    bs = b .+ brng * complexrand(N)
    cs = c .+ crng * complexrand(N)
    zs = z .+ zrng * complexrand(N)

    print("Getting tests... ")
    tests = Vector{NTuple{4, ComplexF64}}()
    tru = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = arb_2f1(ArbComplex.((a,b,c,z), bits = 512)...)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))
    end
    println("done\nTest Count: $(length(tests))")

    # Evaluate each test for accuracy 
    println("\nRunning accuracy tests:")
    print("\tOne Step Taylor: ")
    ta1 = [taylor_2f1(test..., two_step = false)     for test ∈ tests]
    ta1e = clean_error.(ta1, tru)
    println("done")

    print("\tTwo Step Taylor: ")
    ta2 = [taylor_2f1(test..., two_step = true)      for test ∈ tests]
    ta2e = clean_error.(ta2, tru)
    println("done")
    

    # Generate histograms
    fig = Figure()
    ax1 = Axis(fig[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "One Step Taylor")
    ax2 = Axis(fig[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Two Step Taylor")

    bin = 10.0 .^ (-17:1)

    hist!(ax1, ta1e, bins = bin, color = :values, normalization = :probability)
    hist!(ax2, ta2e, bins = bin, color = :values, normalization = :probability)
    
    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))

    resize_to_layout!(fig)

    return fig
end

function reimplot(z, f; kwargs...)
    ref, imf = reim(f)

    fig = Figure()
    axre = Axis3(fig[1,1], title = "Real Part", xlabel = L"Re$(z)$", ylabel = L"Im$(z)$")
    axim = Axis3(fig[1,2], title = "Imaginary Part", xlabel = L"Re$(z)$", ylabel = L"Im$(z)$")

    surface!(axre, reim(z)..., ref; kwargs...)
    surface!(axim, reim(z)..., imf; kwargs...)

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))

    resize_to_layout!(fig)

    return fig
end

function arglines(z, f; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1]; kwargs...)

    x = range(abs.(first(z)), abs.(last(z)), length(z))
    absf = abs.(f)
    arg = angle.(f)
    arg[arg .< 0] .+= π

    plt = lines!(ax, x, absf, 
        colormap = :hsv,
        color = arg
    )

    Colorbar(fig[1,2], plt)

    colsize!(fig.layout, 1, Aspect(1,1))

    resize_to_layout!(fig)

    return fig
end

function surface_error(z, f, title = "Relative Error")
    set_theme!(theme_latexfonts())
    fig = Figure()

    zticks = -15:2:-7
    zticklabels = [latexstring("10^{", val, "}") for val ∈ zticks]

    cticks = -17:2:-5
    cticklabels = [latexstring("10^{", val, "}") for val ∈ cticks]

    ax = Axis3(fig[1,1],
        title = title,
        titlesize = 18,
        xticklabelpad = 0,
        yticklabelpad = 0,
        zticklabelpad = 0,
        xticksvisible = false,
        yticksvisible = false,
        zticksvisible = false,
        xspinesvisible = false,
        yspinesvisible = false,
        zspinesvisible = false,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        perspectiveness = .2,
        elevation = π / 7,
        limits = (nothing, nothing, (-17, -5)),
        zticks = (zticks, zticklabels)
    )

    plt = surface!(ax, reim(z)..., log10.(f),
        colorrange = (-17, -5),
    )

    Colorbar(fig[1,2], plt, ticks = (cticks, cticklabels))

    resize_to_layout!(fig)

    return fig
end

function clean_error(f,t)
    err = abs.((f - t) ./ t)

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-17
        err = 1e-17
    end

    return err
end

function complexrand(N)
    vals = -(1 + im) .+ 2rand(ComplexF64, N)

    return vals
end
