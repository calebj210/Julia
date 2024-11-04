#=
# Tests and graphics for the NSF Grant
#
# Author: Caleb Jacobs
# DLM: November 2, 2024
=#

include("Plotting.jl")
include("pFq.jl")
include("PaperTests/Johansson2F1.jl")
include("PaperTests/Slevinsky2F1.jl")

using BenchmarkTools

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
    z = ComplexGrid(range(-1,3,300), range(-2,2,300))

    println("Running Taylor")
    taylor = average_timings(1, -9/2, -9/4, z, (a,b,c,x) -> taylor_2f1(a, b, c, x; H = .25))
    tp = heatmap(z, taylor; colorscale = log10, colorrange = (1e-7, 1e-2))
    tp.axis.title = L"Taylor Method Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(tp.figure[1,2], tp.plot)
    colsize!(tp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(tp.figure)

    println("Running Levin-Type")
    levin = average_timings(1, -9/2, -9/4, z, weniger_2f1)
    lp = heatmap(z, levin; colorscale = log10, colorrange = (1e-7, 1e-2))
    lp.axis.title = L"Levin-Type Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(lp.figure[1,2], lp.plot)
    colsize!(lp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(lp.figure)

    println("Running Nemo.jl")
    johansson = average_timings(1, -9/2, -9/4, z, johansson_2f1)
    jp = heatmap(z, johansson; colorscale = log10, colorrange = (1e-7, 1e-2))
    jp.axis.title = L"Nemo.jl Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(jp.figure[1,2], jp.plot)
    colsize!(jp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(jp.figure)

    println("Running Mathematica")
    mathematica = average_timings(1, -9/2, -9/4, z, mathematica_2f1)
    mp = heatmap(z, mathematica; colorscale = log10, colorrange = (1e-7, 1e-2))
    mp.axis.title = L"Mathematica Times for ${_2}F_1(1, -9/2; -9/4; z)$"
    Colorbar(mp.figure[1,2], mp.plot)
    colsize!(mp.figure.layout, 1, Aspect(1, 1))
    resize_to_layout!(mp.figure)

    return (; taylor = tp, levin = lp, johansson = jp, mathematica = mp)
end

function average_timings(a, b, c, z, f, N = 5)
    times = zeros(size(z))

    for n ∈ 1 : N
        @time times += [@elapsed f(a, b, c, z) for z ∈ z]
    end
    println("Average time (μs) = ", mean(times) * 1e6)

    times /= N

    return times
end
