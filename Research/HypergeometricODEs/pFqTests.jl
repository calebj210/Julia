#=
# Tests for hypergeometric ODE approach
# 
# Author: Caleb Jacobs
# DLM: September 10, 2024
=#

using CSV, Tables
include("pFq.jl")

# include("Visuals.jl")
using Plots
plotlyjs()

using BenchmarkTools
using MathLink

"Generate CSVs of grid nodes"
function generategrids(path::String, n::Int64, r)
    g = getGrid(n, r)
    CSV.write("Data/" * path, Tables.table([real(g.z) imag(g.z)]), writeheader = false)

    return nothing
end

"Read CSV to complex vector"
function getcomplexvals(path::String)
    reims = CSV.File(path, header = false, delim = ',', types = Float64) |> Tables.matrix

    vals = vec(reims[:, 1] + im * reims[:, 2])

    return vals
end

mathematica_2f1(a, b, c, z) = 
    Complex(weval( W`N[Hypergeometric2F1[a,b,c,z]]`, a = a, b = b, c = c, z = z).args...)
    

"Generate graphics"
function getgraphics(z, f, tru; title = "", dir = -1, exclude = true, logscale = false)
    p1 = complexAbsPlot(z, f - tru, logscale = true, title = title)
    p2 = complexPlot3d(z, f, exclude = exclude, logscale = logscale)
    p3 = complexPlot3d(z, f, T = 2, exclude = exclude, mesh = true, logscale = logscale)
    p4 = complexPlot3d(z, f, T = 3, exclude = exclude, mesh = true, logscale = logscale)
    p5 = complexAbsPlot(z, (f - tru) ./ abs.(tru), logscale = true, title = title)

    camera = attr(
        eye=attr(x = 1.5dir, y = 1.5dir, z = 1.5)
    )
    relayout!(p1, template = "plotly_white")
    relayout!(p2, scene_camera = camera, template = "plotly_white")
    relayout!(p3, scene_camera = camera, template = "plotly_white")
    relayout!(p4, scene_camera = camera, template = "plotly_white")
    relayout!(p5, template = "plotly_white")

    return [p1, p2, p3, p4, p5]
end

function gridtest(a, b, r, n, H = .1, order = 20, taylorN = 100, logscale = false, title = "")
    x = range(-r, r, length = n)
    z = vec(x .+ x' * im)

    CSV.write("Data/grid.csv", Tables.table([real(z) imag(z)]), writeheader = false)
    println("Running grid test with a = ", a, " and b = ", b)
    println("Update Mathematica output and then hit enter to continue.")
    readline()
    tru = getcomplexvals("Data/pfq.csv") 

    F = [fast2f1(a[1], a[2], b[1], z, H = H, order = order, N = taylorN) for z ∈ z]
    
    p = getgraphics(z, F, tru, exclude = true, title = title, logscale = logscale)
    
    return (z, F, tru, p)
end

function convergencetest(a, b, z; order = 20, taylorN = 150, h0 = -2, hf = 0, hN = 100)
    F = Vector{Float64}()
    H = 10 .^ range(h0, hf, length = hN)

    tru = mathematica_2f1(a[1], a[2], b[1], z)

    for h = H
        f = fast2f1(a[1], a[2], b[1], z, H = h, order = order, N = taylorN)
        push!(F, abs(f - tru) / abs(tru))
    end
    
    return (H, F)
end

function timetest(a, b, z; order = 20, N = 150, H = .1)
    t = @benchmark fast2f1($(a[1]), $(a[2]), $(b[1]), $(z), order = $(order), N = $(N), H = $(H))
    val = fast2f1(a[1], a[2], b[1], z, order = order, N = N, H = H)
    tru = mathematica_2f1(a[1], a[2], b[1], z)
    
    med = round(median(t.times) * 1e-3, digits = 2)
    dig = round(Int, -log10(abs(val - tru)))

    return "Median = $(med) μs, Correct Digits = $(dig)"
end

function rungridtests(path::String = ""; N = nothing)
    tests = [
        ([-.9,  1.11], [1.2],    4, 200, .1,  20, 150, false, "2F1(-.9, 1.11; 1.2; z)")

        ([-.9, 5.0-20im], [1.2], 4, 200, .05, 40, 150, false, "2F1(-.9, 5-20im; 1.2; z)")
        ([-.9, 1.11], [50.0im],  4, 200, .05, 40, 150, false, "2F1(-.9, 1.11; 50i; z)")
    ]

    names = ["AbsErr", "AbsArg", "Re", "Im", "RelErr"]

    # Run all tests if not specified
    if isnothing(N)
        N = 1:length(tests)
    end
    
    # Run tests
    p = Vector{Vector{PlotlyJS.SyncPlot}}()
    for (n, test) ∈ pairs(tests)
        if n ∉ N
            continue
        end
        (z, f, tru, figs) = gridtest(test...) 
        push!(p, figs)
        
        # Save figures if wanted
        if path !== ""
            for (f, name) ∈ zip(figs, names)
                savefig(f, path * "Grid$(n)/" * name * ".png", width = 700, height = 700)
            end
        end
    end

    return p
end

function runconvergencetests(path::String = ""; N = nothing, times = false)
    tests = [
        #         a       b        z
        ([ .9, 1.11], [ 1.2], 1.5 + .5im)
        ([ .9, 1.11], [ 1.2], 3cispi(-1/4))
        ([ .9, 1.11], [ 1.2], cispi(1/3)-0.1)

        ([-.9, 1.11], [ 1.2], 1.5 + .5im)
        ([-.9, 1.11], [ 1.2], 3cispi(-1/4))
        ([-.9, 1.11], [ 1.2], cispi(1/3)-0.1)

        ([-0.9, 5.0-20.0im], [1.2], 1.5 - .5im)
        ([-0.9, 5.0-20.0im], [1.2], 3cispi(-1/4))
        ([-0.9, 5.0-20.0im], [1.2], cispi(1/3)-0.1)

        ([-.9, 1.11], [50.0im], 1.5 - .5im)
        ([-.9, 1.11], [50.0im], 3cispi(-1/4))
        ([-.9, 1.11], [50.0im], cispi(1/3)-0.1)
    ]

#     orders = [1,2,4,8,12,16,20,30]
    orders = [30, 40, 60, 100, 150, 200]
    titles = "Test " .* string.(1:length(tests))

    if isnothing(N)
        N = 1:length(tests)
    end

    if !times
        h = []
        plots = Vector{Plots.Plot}()
    end
    for (n, (name, test)) ∈ enumerate(zip(titles, tests))
        if n ∉ N
            continue
        end

        if times
            time = timetest(test..., H = 0.08, order = 30)
            println("Test $(n): $(time)")

            continue
        end

        p = plot(title = name)
        for order ∈ orders
            (h, f) = convergencetest(test..., order = order, hN = 200)
            plot!(h, f, 
                scale  = :log10, 
                xflip  = true, 
                xlabel = "Step Size h",
                ylabel = "Absolute Error",
                ylims  = (1e-17, 1e1),
                yticks = [1e0, 1e-5, 1e-10, 1e-16],
                lw     = 2,
                label  = "T_$(order)")
        end
        plot!(h, h.^orders[1],   ls = :dash, lc = :black, lw = 3, label = "O(h^$(orders[1]))")
        plot!(h, h.^orders[end], ls = :dash, lc = :gray,  lw = 3, label = "O(h^$(orders[end]))")

        push!(plots, p)
    end
    
    if !times
        if path !== ""
            for (n, p) ∈ enumerate(plots)
                savefig(p, path * "Convergence$(n).png")
            end
        end
        return plots
    end
end
