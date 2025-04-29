#=
# Tests and graphics for the NSF Grant
#
# Author: Caleb Jacobs
# DLM: April 27, 2025
=#

include("Plotting.jl")
includet("pFq.jl")
include("PaperTests/Johansson2F1.jl")
include("PaperTests/Slevinsky2F1.jl")
include("../Gregory/HypergeometricGrid.jl")

using ComplexVisuals
using CairoMakie, LaTeXStrings
using BenchmarkTools
using Random
using Suppressor: @suppress_err

# Used Test 1
function levin_errors(a = 1, b = -9/2, c = -9/4)
    z = ComplexGrid(range(-10,10,300), range(-10,10,300))

    print("Getting true solution...")
    tru = johansson_2f1.(a, b, c, z, bits = 106)
    println(" Done")
    print("\tJohansson...")
    johansson = johansson_2f1.(a, b, c, z, bits = 53)
    println(" Done")
    print("\tMathematica...")
    mathematica = mathematica_2f1.(a, b, c, z)
    println(" Done")
    print("\tLevin...")
    levin = weniger_2f1.(a, b, c, z)
    println(" Done")
    print("\tTaylor...")
    # taylor = _2f1.(a, b, c, z)
    taylor = taylor_2f1.(a, b, c, z)
    println(" Done")

    taylor_error = clean_error.(taylor, tru)
    levin_error = clean_error.(levin, tru)
    johansson_error = clean_error.(johansson, tru)
    mathematica_error = clean_error.(mathematica, tru)

    set_theme!(theme_latexfonts())
    ticks = (-15:3:-6, [latexstring("10^{", p, "}") for p ∈ -15:3:-6])

    fig = Figure()
    axt = Axis3(fig[1,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = "Taylor Error",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, (-16,-5)),
    )
    axl = Axis3(fig[1,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = "Levin-Type Error",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, (-16,-5)),
    )
    axj = Axis3(fig[2,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = "Johansson Error",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, (-16,-5)),
    )
    axm = Axis3(fig[2,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Relative Error",
        title = "Mathematica Error",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, (-16,-5)),
    )

    tp = surface!(axt, reim(z)..., log10.(taylor_error),        colorrange = (-16,-6))
    lp = surface!(axl, reim(z)..., log10.(levin_error),         colorrange = (-16,-6))
    jp = surface!(axj, reim(z)..., log10.(johansson_error),     colorrange = (-16,-6))
    mp = surface!(axm, reim(z)..., log10.(mathematica_error),   colorrange = (-16,-6))

    return fig
end

# Used Test 2
function levin_timings(a = 1, b = -9/2, c = -9/4)
    set_theme!(theme_latexfonts())
    z = ComplexGrid(range(-10,10,300), range(-10,10,300))

    println("Running Taylor")
    # taylor = average_time.(a, b, c, z, (a,b,c,x) -> _2f1(a, b, c, x), microseconds = false)
    taylor = average_time.(a, b, c, z, (a,b,c,x) -> taylor_2f1(a, b, c, x), microseconds = false)
    println("Running Levin-Type")
    levin = average_time.(a, b, c, z, weniger_2f1, microseconds = false)
    println("Running Johansson")
    johansson = average_time.(a, b, c, z, (a, b, c, z) -> johansson_2f1(a, b, c, z; bits = 53), microseconds = false)
    println("Running Mathematica")
    mathematica = average_time.(a, b, c, z, mathematica_2f1, microseconds = false)

    mintime = minimum([taylor levin johansson mathematica])
    maxtime = maximum([taylor levin johansson mathematica])
    cbounds = (log10(mintime), log10(maxtime))

    rng = ceil(Int, cbounds[1]):2:floor(Int, cbounds[2])
    ticks = (rng, [latexstring("10^{", p, "}") for p ∈ rng])

    set_theme!(theme_latexfonts())
    fig = Figure()
    axt = Axis3(fig[1,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Time (s)",
        title = "Taylor Time",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, cbounds),
    )
    axl = Axis3(fig[1,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Time (s)",
        title = "Levin-Type Time",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, cbounds),
    )
    axj = Axis3(fig[2,1],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Time (s)",
        title = "Johansson Time",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, cbounds),
    )
    axm = Axis3(fig[2,2],
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = "Time (s)",
        title = "Mathematica Time",
        zticks = ticks,
        xticklabelsize = 12,
        yticklabelsize = 12,
        zticklabelsize = 12,
        zlabelsize = 12,
        xlabeloffset = 25,
        ylabeloffset = 25,
        zlabeloffset = 45,
        limits = (nothing, nothing, cbounds),
    )

    tp = surface!(axt, reim(z)..., log10.(taylor),        colorrange = cbounds)
    lp = surface!(axl, reim(z)..., log10.(levin),         colorrange = cbounds)
    jp = surface!(axj, reim(z)..., log10.(johansson),     colorrange = cbounds)
    mp = surface!(axm, reim(z)..., log10.(mathematica),   colorrange = cbounds)

    # tp = heatmap(z.real, z.imag, taylor; colorscale = log10, colorrange = cbounds)
    # tp.axis.title = latexstring("Taylor Method Times for \${_2}F_1($(a), $(b), $(c); z)\$")
    # tp.axis.xlabel = L"\mathrm{Re}(z)"
    # tp.axis.ylabel = L"\mathrm{Im}(z)"
    # Colorbar(tp.figure[1,2], tp.plot)
    # colsize!(tp.figure.layout, 1, Aspect(1, 1))
    # resize_to_layout!(tp.figure)
    #
    # lp = heatmap(z.real, z.imag, levin; colorscale = log10, colorrange = cbounds)
    # lp.axis.title = latexstring("Levin-Type Times for \${_2}F_1($(a), $(b), $(c); z)\$")
    # lp.axis.xlabel = L"\mathrm{Re}(z)"
    # lp.axis.ylabel = L"\mathrm{Im}(z)"
    # Colorbar(lp.figure[1,2], lp.plot)
    # colsize!(lp.figure.layout, 1, Aspect(1, 1))
    # resize_to_layout!(lp.figure)
    #
    # jp = heatmap(z.real, z.imag, johansson; colorscale = log10, colorrange = cbounds)
    # jp.axis.title = latexstring("Johansson Times for \${_2}F_1($(a), $(b), $(c); z)\$")
    # jp.axis.xlabel = L"\mathrm{Re}(z)"
    # jp.axis.ylabel = L"\mathrm{Im}(z)"
    # Colorbar(jp.figure[1,2], jp.plot)
    # colsize!(jp.figure.layout, 1, Aspect(1, 1))
    # resize_to_layout!(jp.figure)
    #
    # mp = heatmap(z.real, z.imag, mathematica; colorscale = log10, colorrange = cbounds)
    # mp.axis.title = latexstring("Mathematica Times for \${_2}F_1($(a), $(b), $(c); z)\$")
    # mp.axis.xlabel = L"\mathrm{Re}(z)"
    # mp.axis.ylabel = L"\mathrm{Im}(z)"
    # Colorbar(mp.figure[1,2], mp.plot)
    # colsize!(mp.figure.layout, 1, Aspect(1, 1))
    # resize_to_layout!(mp.figure)

    # return (; taylor = tp, levin = lp, johansson = jp, mathematica = mp)
    return fig
end

# Used
function complexrand(N)
    vals = -(1 + im) .+ 2rand(ComplexF64, N)

    return vals
end

# Used Test 3
function random_tests(a = 0, b = 0, c = 0, z = 0; N = 10000, arng = 25, brng = 25, crng = 25, zrng = 1, seed = 997, complextest = false)
    Random.seed!(seed)

    # Setup random tests
    if complextest
        as = a .+ arng * complexrand(N)
        bs = b .+ brng * complexrand(N)
        cs = c .+ crng * complexrand(N)
    else
        as = a .+ arng * (1 .- 2rand(N))
        bs = b .+ brng * (1 .- 2rand(N))
        cs = c .+ crng * (1 .- 2rand(N))
    end
    zs = z .+ zrng * complexrand(N)

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
    println("Test Count: $(length(tests))")

    # Evaluate each test for accuracy 
    print("Running accuracy tests... ")
    ta = [taylor_2f1(test...)                    for test ∈ tests]
    # tr = [_2f1(test...)                          for test ∈ tests]
    le = [weniger_2f1(test...)                   for test ∈ tests]
    jo = [johansson_2f1(test..., bits = 53)      for test ∈ tests]
    ma = @suppress_err [mathematica_2f1(test...) for test ∈ tests]
    print("done\nRunning time tests...\n")

    # Timings of each test
    print("\tTaylor... ")
    tt  = [average_time(test..., taylor_2f1)  for test ∈ tests]
    # print("done\n\tTaylor Transform... ")
    # trt = [average_time(test..., _2f1)        for test ∈ tests]
    print("done\n\tLevin... ")
    lt  = [average_time(test..., weniger_2f1) for test ∈ tests]
    print("done\n\tJohansson... ")
    jt  = [average_time(test..., (a,b,c,z) -> 
                             johansson_2f1(a,b,c,z; bits = 53)) for test ∈ tests]
    print("done\n\tMathematica... ")
    mt = @suppress_err [average_time(test..., mathematica_2f1) for test ∈ tests]
    println("done")

    # Structure data
    taylor      = merge((; avg = round(mean(tt), sigdigits = 2), 
                           med = round(median(tt), sigdigits = 2)), 
                           count_errors(ta, tru))
    # trans       = merge((; avg = round(mean(trt), sigdigits = 2), 
    #                        med = round(median(trt), sigdigits = 2)), 
    #                        count_errors(tr, tru))
    levin       = merge((; avg = round(mean(lt), sigdigits = 2), 
                           med = round(median(lt), sigdigits = 2)), 
                           count_errors(le, tru))
    johansson   = merge((; avg = round(mean(jt), sigdigits = 2), 
                           med = round(median(jt), sigdigits = 2)), 
                           count_errors(jo, tru))
    mathematica = merge((; avg = round(mean(mt), sigdigits = 2), 
                           med = round(median(mt), sigdigits = 2)), 
                           count_errors(ma, tru))

    tae = clean_error.(ta, tru)
    # tre = clean_error.(tr, tru)
    lee = clean_error.(le, tru)
    joe = clean_error.(jo, tru)
    mae = clean_error.(ma, tru)

    # Display result
    print("Taylor:      "); print_error_and_time(taylor)
    # print("Trans:       "); print_error_and_time(trans)
    print("Levin:       "); print_error_and_time(levin)
    print("Johansson:   "); print_error_and_time(johansson)
    print("Mathematica: "); print_error_and_time(mathematica)

    fig = Figure()
    top = GridLayout(fig[1,1])
    bot = GridLayout(fig[2,1])
    axt = Axis(top[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Taylor")
    # axr = Axis(top[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Taylor Transformations")
    axl = Axis(top[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Levin-Type")
    axj = Axis(bot[1,1], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Johansson")
    axm = Axis(bot[1,2], limits = (nothing, (0,1)), xscale = log10, xlabel = "Relative Error", title = "Mathematica")

    bin = 10.0 .^ (-16:2:2)

    hist!(axt, tae, bins = bin, color = :values, normalization = :probability)
    # hist!(axr, tre, bins = bin, color = :values, normalization = :probability)
    hist!(axl, lee, bins = bin, color = :values, normalization = :probability)
    hist!(axj, joe, bins = bin, color = :values, normalization = :probability)
    hist!(axm, mae, bins = bin, color = :values, normalization = :probability)
    
    rowsize!(fig.layout, 1, Auto(1))
    rowsize!(fig.layout, 2, Auto(1))

    return fig
end

# Used
function count_errors(vals, true_values; good = 1e-14, fair = 1e-9, poor = 1e-1)    
    goods = fairs = poors = fails = 0

    errs = abs.((vals - true_values) ./ true_values)

    for (err, tru) ∈ zip(errs, true_values)
        if isnan(tru) || isinf(tru)
            continue
        end
        if err <= good
            goods += 1
        elseif err <= fair
            fairs += 1
        elseif err <= poor
            poors += 1
        else
            fails += 1
        end
    end

    return (; goods, fairs, poors, fails)
end

# Used
function print_error_and_time(r)
    print(r.avg, " Average, ", r.med, " Median, ", r.goods, " Good, ", r.fairs, " Fair, ", r.poors, " Poor, ", r.fails, " Fails.\n")
end

# Used
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

# Used
function mathematica_grid_error()
    z = complex_square_grid(2.5,300)
    tru = johansson_2f1.(1.99,.9,2.9,z,bits = 256)
    mathematica = mathematica_2f1.(1.99,.9,2.9,z)

    err = abs.((mathematica - tru) ./ tru)
    err[iszero.(err)] .= 1e-17

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(fig[1,1], 
              xlabel = L"\mathrm{Re}(z)",
              ylabel = L"\mathrm{Im}(z)",
              title = L"Mathematica Error in ${_2}F_1(1.99, 0.9; 2.9; z)$"
             )

    plt = heatmap!(ax, z, err,
             colorscale = log10
            )
    Colorbar(fig[1,2], plt)
    colsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    return fig
end

# Used
function clean_error(f,t)
    err = abs.((f - t) ./ t)

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-16
        err = 1e-16
    end

    return err
end

function taylor_terms(a, b, c, z, p, step_max = Inf, init_max = exp(-1))
    if real(z) > 1
        z0 = imag(z) > 0 ? im * init_max : -im * init_max
    else
        z0 = sign(z) * init_max
    end
    dir = sign(z - z0)

    h_opt = abs(z0 - 1) / exp(2)
    # h_opt = abs(z0 - 1)
    h_end = abs(z0 - z)
    # h_ord = sqrt(abs(2zn * (1 - z0) / (a * b)))                    # sqrt(C / A)
    # h_ord = abs(2zn * (1 - z0) / (c - (1 + a + b) * z0))                  # C
    # h_ord = abs(2zn * (1 - z0) * (c - (1 + a + b) * z0) / (a * b))        # C B / A
    h_ord = Inf

    h = dir * min(h_opt, h_end, h_ord, step_max)        # Step size based on Jorba and Zou 2005

    coeffs = taylor_coefficients(a, b, c, z0, p)

    terms = coeffs .* h.^(0:(p - 1))

    fig = Figure()
    axt = Axis(fig[1,1], yscale = log10)
    axs = Axis(fig[2,1], yscale = log10)
    axe = Axis(fig[3,1], yscale = log10)
    lines!(axt, abs.(coeffs), label = L"|a_n|")
    lines!(axs, abs.(cumsum(terms)), label = L"|S_n|")
    lines!(axe, abs.(terms), label = L"|a_n h^n|")
    lines!(axe, eps.(abs.(cumsum(terms))), label = L"\varepsilon(|S_n|)")
    Legend(fig[1,2], axt)
    Legend(fig[2,2], axs)
    Legend(fig[3,2], axe)
    resize_to_layout!(fig)

    init = maclaurin_2f1(a,b,c,z0)[1]
    tru_init = johansson_2f1(a,b,c,z0,bits = 1024)
    err = abs.((init - tru_init) ./ tru_init)

    print("Test (a,b,c,z):\n\ta = $(a)\n\tb = $(b)\n\tc = $(c)\n\tz = $(z)\n")
    print("Initialization error = $(err)\n")

    return (fig, terms)
end

function slevinsky_comparison()
    z = ComplexGrid(range(-10,10,300), range(-10,10,300))

    println("Johansson Time")
    johansson = johansson_2f1.(1, -9/2, -9/4, z)
    @btime johansson = johansson_2f1.(1, -9/2, -9/4, $(z))

    println("Levin-Type Time")
    levin = weniger_2f1.(1, -9/2, -9/4, z)
    @btime levin = weniger_2f1.(1, -9/2, -9/4, $(z))

    println("Taylor Time")
    taylor = taylor_2f1.(1, -9/2, -9/4, z; H = .25)
    @btime taylor = taylor_2f1.(1, -9/2, -9/4, $(z); H = .25)

    levin_plots = grid_error_plot(z, levin, johansson)
    levin_plots.axis.title = L"Levin-Type for ${_2}F_1(1, -9/2; -9/4; z)$"
    println("Levin Max-Relative-Error = ", maximum(abs.((levin - johansson) ./ johansson)))

    taylor_plots = grid_error_plot(z, taylor, johansson)
    taylor_plots.axis.title = L"Taylor Method for ${_2}F_1(1, -9/2; -9/4; z)$"
    println("Taylor Max-Relative-Error = ", maximum(abs.((taylor - johansson) ./ johansson)))

    return (;taylor = taylor_plots, levin = levin_plots)
end

function average_grid_timings()
    Ns = 10:10:300
    (a, b, c) = (1, -9/2, -9/4)

    taylor_times        = zeros(length(Ns))
    taylor_means        = zeros(length(Ns))
    taylor_medians      = zeros(length(Ns))
    taylor_l2           = zeros(length(Ns))

    levin_times         = zeros(length(Ns))
    levin_means         = zeros(length(Ns))
    levin_medians       = zeros(length(Ns))
    levin_l2            = zeros(length(Ns))

    johansson_times     = zeros(length(Ns))
    johansson_means     = zeros(length(Ns))
    johansson_medians   = zeros(length(Ns))
    johansson_l2        = zeros(length(Ns))

    mathematica_times   = zeros(length(Ns))
    mathematica_means   = zeros(length(Ns))
    mathematica_medians = zeros(length(Ns))
    mathematica_l2      = zeros(length(Ns))

    grid_spacings       = zeros(length(Ns))

    for (i, N) ∈ enumerate(Ns)
        z = ComplexGrid(range(-1,3,N), range(-2,2,N))

        println("Test ", i, " of ", length(Ns), ":")
        tru = johansson_2f1.(a, b, c, z, bits = 128)

        ft = zeros(ComplexF64, N, N)
        fl = zeros(ComplexF64, N, N)
        fj = zeros(ComplexF64, N, N)
        fm = zeros(ComplexF64, N, N)

        print("\tTaylor:      ")
        taylor_times[i] = @elapsed ft = taylor_2f1.(a, b, c, z)
        println(round(taylor_times[i], sigdigits = 3), "s") 

        print("\tLevin-Type:  ")
        @suppress_err levin_times[i] = @elapsed fl = weniger_2f1.(a, b, c, z)
        println(round(levin_times[i], sigdigits = 3), "s") 

        print("\tJohansson:   ")
        johansson_times[i] = @elapsed fj = johansson_2f1.(a, b, c, z, bits = 53)
        println(round(johansson_times[i], sigdigits = 3), "s") 

        print("\tMathematica: ")
        mathematica_times[i] = @elapsed fm = mathematica_2f1.(a, b, c, z)
        println(round(mathematica_times[i], sigdigits = 3), "s") 

        relerr = abs.((ft - tru) ./ tru)
        taylor_means[i]   = mean(relerr)
        taylor_medians[i] = median(relerr)
        taylor_l2[i] = norm(ft - tru) / norm(tru)

        relerr = abs.((fl - tru) ./ tru)
        levin_means[i]   = mean(relerr)
        levin_medians[i] = median(relerr)
        levin_l2[i] = norm(fl - tru) / norm(tru)

        relerr = abs.((fj - tru) ./ tru)
        johansson_means[i]   = mean(relerr)
        johansson_medians[i] = median(relerr)
        johansson_l2[i] = norm(fj - tru) / norm(tru)

        relerr = abs.((fm - tru) ./ tru)
        mathematica_means[i]   = mean(relerr)
        mathematica_medians[i] = median(relerr)
        mathematica_l2[i] = norm(fm - tru) / norm(tru)

        grid_spacings[i] = abs(z[1] - z[2])
    end

    times = (; taylor_times, levin_times, johansson_times, mathematica_times)
    errs  = (; grid_spacings,
               taylor_means,      taylor_medians,      taylor_l2,
               levin_means,       levin_medians,       levin_l2,
               johansson_means,   johansson_medians,   johansson_l2,
               mathematica_means, mathematica_medians, mathematica_l2)

    return (; times, errs)
end

function plot_average_grid_timings(times, errs)
    set_theme!(theme_latexfonts())
    fig = Figure()
    
    ax1 = Axis(fig[1,1], xscale = log10, yscale = log10, xreversed = true)
    ax1.ylabel = "Time (seconds)"
    ax1.title = L"Time Complexity and Error for ${_2}F_1$"
    ax1.xminorticksvisible = true
    ax1.xminorticks = IntervalsBetween(5)
    hidexdecorations!(ax1, minorgrid = false, grid = false)

    ax2 = Axis(fig[2,1], xscale = log10, yscale = log10, xreversed = true)
    ax2.xlabel = "Grid Spacing"
    ax2.ylabel = "Relative Error"

    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 10 
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace

    linkxaxes!(ax1, ax2)

    lines!(ax1, errs.grid_spacings, times.taylor_times;      color = :orange, label = "Taylor")
    lines!(ax1, errs.grid_spacings, times.levin_times;       color = :blue,   label = "Levin-Type")
    lines!(ax1, errs.grid_spacings, times.johansson_times;   color = :green,  label = "Johansson")
    lines!(ax1, errs.grid_spacings, times.mathematica_times; color = :purple, label = "Mathematica")

    lines!(ax2, errs.grid_spacings, errs.taylor_means;      color = :orange, linestyle = :solid, label = "Mean")
    lines!(ax2, errs.grid_spacings, errs.levin_means;       color = :blue,   linestyle = :solid, label = "Mean")
    lines!(ax2, errs.grid_spacings, errs.johansson_means;   color = :green,  linestyle = :solid, label = "Mean")
    lines!(ax2, errs.grid_spacings, errs.mathematica_means; color = :purple, linestyle = :solid, label = "Mean")

    lines!(ax2, errs.grid_spacings, errs.taylor_medians;      color = :orange, linestyle = :dash, label = "Median")
    lines!(ax2, errs.grid_spacings, errs.levin_medians;       color = :blue,   linestyle = :dash, label = "Median")
    lines!(ax2, errs.grid_spacings, errs.johansson_medians;   color = :green,  linestyle = :dash, label = "Median")
    lines!(ax2, errs.grid_spacings, errs.mathematica_medians; color = :purple, linestyle = :dash, label = "Median")

    lines!(ax2, errs.grid_spacings, errs.taylor_l2;      color = :orange, linestyle = :dot, label = L"$L^2$ Norm")
    lines!(ax2, errs.grid_spacings, errs.levin_l2;       color = :blue,   linestyle = :dot, label = L"$L^2$ Norm")
    lines!(ax2, errs.grid_spacings, errs.johansson_l2;   color = :green,  linestyle = :dot, label = L"$L^2$ Norm")
    lines!(ax2, errs.grid_spacings, errs.mathematica_l2; color = :purple, linestyle = :dot, label = L"$L^2$ Norm")

    fig[1,2] = Legend(fig, ax1, "Method", framevisible = false)

    means   = LineElement(color = :black, linestyle = :solid)
    medians = LineElement(color = :black, linestyle = :dash)
    l2s     = LineElement(color = :black, linestyle = :dot)

    Legend(fig[2,2], [means, medians, l2s], ["Mean", "Median", L"$L^2$ Norm"], "Error Type", framevisible = false)

    return fig
end

function run_grid_complexity_tests(a, b, c; r = 1.99, Ns = 81:22:301)
    fig = Figure()
    ax = Axis(fig[1,1]; xscale = log10, yscale = log10)

    println("Taylor running")
    grid_complexity_test!(ax, a, b, c, taylor_2f1; r = r, Ns = Ns, label = "Taylor")
    println("Levin running")
    @suppress_err grid_complexity_test!(ax, a, b, c, weniger_2f1; r = r, Ns = Ns, label = "Levin")
    println("Johansson running")
    grid_complexity_test!(ax, a, b, c, (a,b,c,z) -> johansson_2f1(a,b,c,z; bits = 53); r = r, Ns = Ns, label = "Johansson")
    println("Mathematica running")
    grid_complexity_test!(ax, a, b, c, mathematica_2f1; r = r, Ns = Ns, label = "Mathematica")
    # println("End-corrected trapezoidal rule running")
    # grid_complexity_test!(ax, a, b, c; r = r, Ns = Ns, label = "Trapezoid")

    fig[1,2] = Legend(fig, ax, "Methods", framevisible = false)

    return Makie.FigureAxis(fig, ax)
end

function grid_complexity_test!(ax, a, b, c, f; r = 2, Ns = 81:22:301, kwargs...)
    times = zeros(length(Ns))

    for (i, N) ∈ enumerate(Ns)
        z = complex_square_grid(r, N)

        times[i] = @elapsed f.(a, b, c, z)
    end

    lines!(ax, Ns, times; kwargs...)
end

function grid_complexity_test!(ax, a, b, c; r = 2, Ns = 81:22:301, kwargs...)
    times = zeros(length(Ns))

    for (i, N) ∈ enumerate(Ns)
        times[i] = @elapsed pFq([a, b], [c]; 
                                grid_radius = r,
                                grid_points = (N - 1) ÷ 2,
                                padding_layers = 5,
                                taylor_radius = .5,
                                modify_z1 = true,
                                correction_radius = .5,
                                inner_radius = .6,
                                outer_radius = .8,
                                z1_expansion_order = 70)
    end

    lines!(ax, Ns, times; kwargs...)
end

function trapezoidal_rule_test(a, b; r = 2.49, Ns = 41:4:101)
    trapezoid_times   = zeros(length(Ns))
    trapezoid_means   = zeros(length(Ns))
    trapezoid_medians = zeros(length(Ns))
    trapezoid_l2      = zeros(length(Ns))

    levin_times       = zeros(length(Ns))
    levin_means       = zeros(length(Ns))
    levin_medians     = zeros(length(Ns))
    levin_l2          = zeros(length(Ns))

    mathematica_times = zeros(length(Ns))

    grid_spacings     = zeros(length(Ns))

    for (i, N) ∈ enumerate(Ns)
        print("Trapezoidal Test ", i, " of ", length(Ns), ": ")
        trapezoid_times[i] = @elapsed (z,ft,_) = pFq(a, b; 
                                          grid_radius = r,
                                          grid_points = N,
                                          padding_layers = 5,
                                          taylor_radius = .5,
                                          modify_z1 = true,
                                          correction_radius = .5,
                                          inner_radius = .6,
                                          outer_radius = .8,
                                          z1_expansion_order = 120)
        println(round(trapezoid_times[i], sigdigits = 3), "s") 
        
        fl = copy(ft)
        fm = copy(ft)

        print("Levin-Type  Test ", i, " of ", length(Ns), ": ")
        @suppress_err levin_times[i] = @elapsed fl = [weniger_pfq((a...,), (b...,), z) for z ∈ z]
        println(round(levin_times[i], sigdigits = 3), "s") 

        print("Mathematica Test ", i, " of ", length(Ns), ": ")
        mathematica_times[i] = @elapsed fm = [mathematica_pfq(a, b, z) for z ∈ z]
        println(round(mathematica_times[i], sigdigits = 3), "s") 

        relerr = abs.((ft - fm) ./ fm)
        trapezoid_means[i]   = mean(relerr)
        trapezoid_medians[i] = median(relerr)
        trapezoid_l2[i] = norm(ft - fm) / norm(fm)

        relerr = abs.((fl - fm) ./ fm)
        levin_means[i]   = mean(relerr)
        levin_medians[i] = median(relerr)
        levin_l2[i] = norm(fl - fm) / norm(fm)

        grid_spacings[i] = abs(z[1] - z[2])
    end

    fig = Figure()
    ax = Axis(fig[1,1], xscale = log10, yscale = log10, xreversed = true)
    ax.xlabel = "Grid Spacing"
    ax.ylabel = "Time (seconds)"

    lines!(ax, grid_spacings, trapezoid_times;   label = "Trapezdoidal")
    lines!(ax, grid_spacings, levin_times;       label = "Levin-Type")
    lines!(ax, grid_spacings, mathematica_times; label = "Mathematica")

    fig[1,2] = Legend(fig, ax, "Method", framevisible = false)

    return (fig, ax, 
            (; trapezoid_times, levin_times, mathematica_times), 
            (;grid_spacings, 
             trapezoid_means, trapezoid_medians, trapezoid_l2,
             levin_means,     levin_medians,     levin_l2)
           )
end

function make_time_error_figure(times, errs)
    set_theme!(theme_latexfonts())
    fig = Figure()
    ax1 = Axis(fig[1,1], xscale = log10, yscale = log10, xreversed = true)
    ax1.xminorticksvisible = true
    ax1.xminorticks = IntervalsBetween(5)
    hidexdecorations!(ax1, minorgrid = false, grid = false)
    ax1.ylabel = "Time (seconds)"

    ax2 = Axis(fig[2,1], xscale = log10, yscale = log10, xreversed = true)
    ax2.xminorticksvisible = true
    ax2.xminorticks = IntervalsBetween(5)
    ax2.xlabel = "Grid Spacing"
    ax2.ylabel = "Relative Error"

    linkxaxes!(ax1, ax2)

    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 10 
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace

    lines!(ax1, errs.grid_spacings, times.trapezoid_times;   color = :orange, label = "Trapezdoidal")
    lines!(ax1, errs.grid_spacings, times.levin_times;       color = :blue,   label = "Levin-Type")
    lines!(ax1, errs.grid_spacings, times.mathematica_times; color = :green,  label = "Mathematica")

    lines!(ax2, errs.grid_spacings, errs.trapezoid_means;   color = :orange, linestyle = :solid, label = "Mean")
    lines!(ax2, errs.grid_spacings, errs.levin_means;       color = :blue,   linestyle = :solid, label = "Mean")
    lines!(ax2, errs.grid_spacings, errs.trapezoid_medians; color = :orange, linestyle = :dash,  label = "Median")
    lines!(ax2, errs.grid_spacings, errs.levin_medians;     color = :blue,   linestyle = :dash,  label = "Median")
    lines!(ax2, errs.grid_spacings, errs.trapezoid_l2;      color = :orange, linestyle = :dot,   label = L"$L^2$ Norm")
    lines!(ax2, errs.grid_spacings, errs.levin_l2;          color = :blue,   linestyle = :dot,   label = L"$L^2$ Norm")

    fig[1,2] = Legend(fig, ax1, "Method", framevisible = false)

    means   = LineElement(color = :black, linestyle = :solid)
    medians = LineElement(color = :black, linestyle = :dash)
    l2s     = LineElement(color = :black, linestyle = :dot)

    Legend(fig[2,2], [means, medians, l2s], ["Mean", "Median", L"$L^2$ Norm"], "Error Type", framevisible = false)

    ax1.title = L"Time Complexity and Error for ${_5}F_4$"

    return fig
end

function crespo_tests_table_1()
    tests = [
        (-.1,.2,.3,.5)
        (-.1,.2,.3,1.5)
        (-.1,.2,.3,100)
        (2+8im,3-5im,sqrt(2) - im*π,.25)
        (2+8im,3-5im,sqrt(2) - im*π,.75)
        (2+8im,3-5im,sqrt(2) - im*π,-10)
        (2+200im,5-100im,10+500im,.8)
        (2.25,3.75,-.5,-1)
        (2+200im,5,10,0.6)
        (1,3,7,.25)
        (1,3,7,.75)
        (1,3,7,-3)
       ]

    johansson = [johansson_2f1(test...) for test ∈ tests]
    taylor = [taylor_2f1(test...) for test ∈ tests]

    relative_error = round.(abs.((taylor - johansson) ./ johansson), sigdigits = 2)

    return relative_error
end

function crespo_tests_table_2()
    tests = [
        (.1,.2,-.3,-.5+.5im)
        (.1,.2,-.3,1+.5im)
        (.1,.2,-.3,5+5im)
        (4,1.1,2,cispi(1/3))
        (4,1.1,2,1+5im)
        (4,1.1,2,-5+5im)
        (2/3,1,4/3,cispi(1/3))
        (2/3,1,4/3,2im)
        (2/3,1,4/3,1+im)
        (2/3,1,4/3,100im)
       ]

    johansson = [johansson_2f1(test...) for test ∈ tests]
    taylor = [taylor_2f1(test...) for test ∈ tests]

    relative_error = round.(abs.((taylor - johansson) ./ johansson), sigdigits = 2)

    return relative_error
end
