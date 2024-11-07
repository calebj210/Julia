#=
# Tests and graphics for the NSF Grant
#
# Author: Caleb Jacobs
# DLM: November 6, 2024
=#

include("Plotting.jl")
include("pFq.jl")
include("PaperTests/Johansson2F1.jl")
include("PaperTests/Slevinsky2F1.jl")
include("../Gregory/HypergeometricGrid.jl")

using CairoMakie
using BenchmarkTools
using Suppressor: @suppress_err

function levin_error()
    z = ComplexGrid(range(-1,3,300), range(-2,2,300))

    tru = johansson_2f1.(1, -9/2, -9/4, z; bits = 106)
    levin = weniger_2f1.(1, -9/2, -9/4, z)
    taylor = taylor_2f1.(1, -9/2, -9/4, z)

    levin_error = abs.((levin - tru) ./ tru)
    levin_error[iszero.(levin_error)] .= 1e-17
    taylor_error = abs.((taylor - tru) ./ tru)
    taylor_error[iszero.(taylor_error)] .= 1e-17

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel = "Re(z)", ylabel = "Im(z)", title = L"Taylor Method Error for ${_2}F_1(1, -9/4; -9/2; z)$")
    ax2 = Axis(fig[1,2], xlabel = "Re(z)", title = L"Taylor Method Error for ${_2}F_1(1, -9/4; -9/2; z)$")

    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)
    
    plt1 = heatmap!(ax1, z.real, z.imag, taylor_error; colorrange = (1e-17, 1e-5), colorscale = log10)
    plt2 = heatmap!(ax2, z.real, z.imag, levin_error;  colorrange = (1e-17, 1e-5), colorscale = log10)

    Colorbar(fig[1,3], plt1)
    colsize!(fig.layout, 1, Aspect(1, 1))
    colsize!(fig.layout, 2, Aspect(1, 1))

    resize_to_layout!(fig)

    return fig
end

function slevinsky_comparison()
    z = ComplexGrid(range(-1,3,300), range(-2,2,300))

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

function grid_timings()
    set_theme!(theme_latexfonts())
    z = ComplexGrid(range(-1,3,300), range(-2,2,300))

    println("Running Taylor")
    taylor = average_time.(1, -9/2, -9/4, z, (a,b,c,x) -> taylor_2f1(a, b, c, x))
    tp = heatmap(z.real, z.imag, taylor; colorscale = log10, colorrange = (1e-7, 1e-2))
    tp.axis.title = L"Taylor Method Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(tp.figure[1,2], tp.plot)
    colsize!(tp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(tp.figure)

    println("Running Levin-Type")
    levin = average_time.(1, -9/2, -9/4, z, weniger_2f1)
    lp = heatmap(z.real, z.imag, levin; colorscale = log10, colorrange = (1e-7, 1e-2))
    lp.axis.title = L"Levin-Type Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(lp.figure[1,2], lp.plot)
    colsize!(lp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(lp.figure)

    println("Running Johansson")
    johansson = average_time.(1, -9/2, -9/4, z, (a, b, c, z) -> johansson_2f1(a, b, c, z; bits = 106))
    jp = heatmap(z.real, z.imag, johansson; colorscale = log10, colorrange = (1e-7, 1e-2))
    jp.axis.title = L"Johansson Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(jp.figure[1,2], jp.plot)
    colsize!(jp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(jp.figure)

    println("Running Mathematica")
    mathematica = average_time.(1, -9/2, -9/4, z, mathematica_2f1)
    mp = heatmap(z.real, z.imag, mathematica; colorscale = log10, colorrange = (1e-7, 1e-2))
    mp.axis.title = L"Mathematica Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(mp.figure[1,2], mp.plot)
    colsize!(mp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(mp.figure)

    return (; taylor = tp, levin = lp, johansson = jp, mathematica = mp)
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

function average_time(a, b, c, z, f, N = 5)
    time = 0.0

    for n ∈ 1 : N
        time += @elapsed f(a, b, c, z)
    end

    time /= N

    time = round(Int, time * 1e6) # Convert to microseconds

    return time
end

function count_errors(vals, true_values; good = 1e-14, fair = 1e-9, poor = 1e-1)
    goods = fairs = poors = fails = 0

    errs = abs.((vals - true_values) ./ true_values)

    for err ∈ errs
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

function random_tests(a = 1, b = 2, c = 3, z = .5 + .5im; N = 10001, rng = 10, seed = 997)
    Random.seed!(seed)

    # Setup random tests
    as = a .+ 2rng * (rand(ComplexF64, N) .- .5(1 + im))
    bs = b .+ 2rng * (rand(ComplexF64, N) .- .5(1 + im))
    cs = c .+ 2rng * (rand(ComplexF64, N) .- .5(1 + im))
    zs = z .+ 2rng * (rand(ComplexF64, N) .- .5(1 + im))
    collected_tests = [(a, b, c, z) for (a,b,c) ∈ zip(as, bs, cs)]

    # Evaluate each test for accuracy 
    true_vals        = [johansson_2f1(test...)   for test ∈ collected_tests]
    taylor_vals      = [taylor_2f1(test...)      for test ∈ collected_tests]
    levin_vals       = [weniger_2f1(test...)     for test ∈ collected_tests]
    mathematica_vals = [mathematica_2f1(test...) for test ∈ collected_tests]
    johansson_vals   = [johansson_2f1(test...,
                                     bits = 53)  for test ∈ collected_tests]

    # Timings of each test
    true_times        = [average_time(test..., johansson_2f1)   for test ∈ collected_tests]
    taylor_times      = [average_time(test..., taylor_2f1)      for test ∈ collected_tests]
    levin_times       = [average_time(test..., weniger_2f1)     for test ∈ collected_tests]
    mathematica_times = [average_time(test..., mathematica_2f1) for test ∈ collected_tests]
    johansson_times   = [average_time(test..., (a,b,c,z) -> 
                             johansson_2f1(a,b,c,z; bits = 53)) for test ∈ collected_tests]

    # Structure data
    true_jo     = merge((; avg = round(mean(true_times), sigdigits = 2), 
                           med = round(median(true_times), sigdigits = 2)), 
                           count_errors(true_vals, true_vals))
    taylor      = merge((; avg = round(mean(taylor_times), sigdigits = 2), 
                           med = round(median(taylor_times), sigdigits = 2)), 
                           count_errors(taylor_vals, true_vals))
    levin       = merge((; avg = round(mean(levin_times), sigdigits = 2), 
                           med = round(median(levin_times), sigdigits = 2)), 
                           count_errors(levin_vals, true_vals))
    johansson   = merge((; avg = round(mean(johansson_times), sigdigits = 2), 
                           med = round(median(johansson_times), sigdigits = 2)), 
                           count_errors(johansson_vals, true_vals))
    mathematica = merge((; avg = round(mean(mathematica_times), sigdigits = 2), 
                           med = round(median(mathematica_times), sigdigits = 2)), 
                           count_errors(mathematica_vals, true_vals))

    # Display result
    print("True:        "); print_error_and_time(true_jo)
    print("Taylor:      "); print_error_and_time(taylor)
    print("Levin:       "); print_error_and_time(levin)
    print("Johansson:   "); print_error_and_time(johansson)
    print("Mathematica: "); print_error_and_time(mathematica)

    return (;taylor, levin, mathematica, johansson)
end

function print_error_and_time(r)
    print(r.avg, " Average, ", r.med, " Median, ", r.goods, " Good, ", r.fairs, " Fair, ", r.poors, " Poor, ", r.fails, " Fails.\n")
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
