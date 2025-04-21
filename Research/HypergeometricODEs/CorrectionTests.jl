include("pFq.jl")
using ComplexVisuals, CairoMakie, LaTeXStrings
import Random.seed!

function correction_test()
    a, b, c = (-19.139483726210333, -14.134535062801518, -24.30112282763963)

    fig = Figure()
    ax1 = Axis3(fig[1,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"\log_{10}|z|",
        title = L"f_1(z)",
        elevation = .2π,
        viewmode = :fit,
        )
    ax2 = Axis3(fig[1,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"\log_{10}|z|",
        title = L"f_2(z)",
        elevation = .2π,
        viewmode = :fit,
    )
    
    z = complex_square_grid(2, 300)
    f1 = johansson_2f1.(a, b, c, z)
    f2 = z.^(1 - c) .* johansson_2f1.(a - c + 1, b - c + 1, 2 - c, z)

    complexsurface!(ax1, z, Makie.pseudolog10.(f1))
    complexsurface!(ax2, z, Makie.pseudolog10.(f2))

    rowsize!(fig.layout, 1, Aspect(1,1))

    resize_to_layout!(fig)

    return fig
end

function alternate_error(z = nothing, f1 = nothing, f2 = nothing)
    a, b, c = (-19.139483726210333, -14.134535062801518, -24.30112282763963)

    ticks = (-15:3:0, [latexstring("10^{", p, "}") for p ∈ -15:3:0])

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax1 = Axis3(fig[1,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = L"Error in $f_1(z)$",
        titlegap = -12,
        elevation = .2π,
        viewmode = :fit,
        limits = (nothing,nothing,(-16,1)),
        zticks = ticks,
    )
    ax2 = Axis3(fig[1,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = L"Error in $f_2(z)$",
        titlegap = -12,
        elevation = .2π,
        viewmode = :fit,
        limits = (nothing,nothing,(-16,1)),
        zticks = ticks,
    )

    if isnothing(z) || isnothing(f1) || isnothing(f2)
        z = complex_square_grid(2, 300)

        f1 = johansson_2f1.(a, b, c, z)
        f2 = z.^(1 - c) .* johansson_2f1.(a - c + 1, b - c + 1, 2 - c, z)
    end

    t1 = taylor_2f1.(a, b, c, z)
    t2 = z.^(1 - c) .* taylor_2f1.(a - c + 1, b - c + 1, 2 - c, z)

    e1 = clean_error.(t1, f1)
    e2 = clean_error.(t2, f2)

    surface!(ax1, reim(z)..., log10.(e1),
             colorrange = (-16, 1),
    )
    surface!(ax2, reim(z)..., log10.(e2),
             colorrange = (-16, 1),
    )

    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))

    Colorbar(fig[1,3], limits = (-16, 1), ticks = ticks)

    resize_to_layout!(fig)

    return ((z, f1, f2), fig)
end

function random_failed_tests(a = 0, b = 0, c = 0, z = 0; N = 10000, arng = 25, brng = 25, crng = 25, zrng = 1, seed = 997)
    seed!(seed)

    # Setup random tests
    as = a .+ arng * (1 .- 2rand(N))
    bs = b .+ brng * (1 .- 2rand(N))
    cs = c .+ crng * (1 .- 2rand(N))
    # zs = z .+ zrng * (1 .- rand(N))
    # as = a .+ arng * complexrand(N)
    # bs = b .+ brng * complexrand(N)
    # cs = c .+ crng * complexrand(N)
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
    print("\tRegular Taylor: ")
    ta = [taylor_2f1(test..., backward = false) for test ∈ tests]
    tae = clean_error.(ta, tru)
    print("done\n\tBackwards Taylor: ")
    tc = [taylor_2f1(test..., backward = true)  for test ∈ tests]
    tce = clean_error.(tc, tru)
    println("done")
    
    print("Collecting tests: ")
    taf = tests[isone.(tae)]
    tcf = tests[isone.(tce)]
    println("done")

    # Generate histograms
    fig = Figure()
    axt = Axis(fig[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Regular Taylor")
    axc = Axis(fig[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Backwards Corrected Taylor")

    bin = 10.0 .^ (-16:2:2)

    hist!(axt, tae, bins = bin, color = :values, normalization = :probability)
    hist!(axc, tce, bins = bin, color = :values, normalization = :probability)
    
    colsize!(fig.layout, 1, Aspect(1,1))
    colsize!(fig.layout, 2, Aspect(1,1))

    resize_to_layout!(fig)

    return (; fig, taf, tcf)
end

function clean_error(f,t)
    err = abs.((f - t) ./ t)

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
