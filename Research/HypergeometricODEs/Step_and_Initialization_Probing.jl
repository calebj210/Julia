#= 
# Tests for inspecting failures of the Taylor method
#
# Author: Caleb Jacobs
# DLM: March 4, 2025
=#

include("pFq.jl")
include("PaperTests/Slevinsky2F1.jl")

using CairoMakie, ComplexVisuals
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
    val = taylor_2f1.(a, b, c, z)
    println("done!")

    errs = clean_error.(val, tru)

    set_theme!(complex_theme)
    fig = Figure()

    ax1 = Axis(fig[1,1], title = "Error")
    ax2 = Axis(fig[1,3], title = "Taylor Phase")
    ax3 = Axis(fig[1,4], title = "True Phase")

    plt1 = heatmap!(ax1, z, errs, colorscale = log10, colorrange = (1e-17,1))
    plt2 = phase!(ax2, z, val)
    plt3 = phase!(ax3, z, tru)

    Colorbar(fig[1,2], plt1)

    for n ∈ [1,3,4]
        colsize!(fig.layout, n, Aspect(1, 1))
    end
    resize_to_layout!(fig)

    return (fig, z, tru)
end

function random_failed_tests(a = 0, b = 0, c = 0, z = 0; N = 10000, arng = 30, brng = 30, crng = 30, zrng = 1, seed = 997)
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
    ta = [taylor_2f1(test...)               for test ∈ tests]
    tae = clean_error.(ta, tru)
    taf = tests[isone.(tae)]
    print("done\n\tTransforms: ")
    tr = [_2f1(test...)                     for test ∈ tests]
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
