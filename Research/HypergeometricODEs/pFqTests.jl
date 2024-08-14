#=
# Tests for hypergeometric ODE approach
# 
# Author: Caleb Jacobs
# DLM: August 14, 2024
=#

using CSV, Tables
# include("Visuals.jl")
using Plots
plotlyjs()
include("pFq.jl")

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

function gridtest(a, b, r, n, h = .1, order = 20, taylorN = 100, logscale = false, title = "")
    x = range(-r, r, length = n)
    z = vec(x .+ x' * im)

    CSV.write("Data/grid.csv", Tables.table([real(z) imag(z)]), writeheader = false)
    println("Running grid test with a = ", a, " and b = ", b)
    println("Update Mathematica output and then hit enter to continue.")
    readline()
    tru = getcomplexvals("Data/pfq.csv") 

    F = [F21(a[1], a[2], b[1], z, h = h, order = order, taylorN = taylorN) for z ∈ z]
    
    p = getgraphics(z, F, tru, exclude = true, title = title, logscale = logscale)
    
    return (z, F, tru, p)
end

function convergencetest(a, b, z, tru; z0 = 0, order = 20, taylorN = 100, h0 = -2, hf = 0, hN = 100)
    F = Vector{Float64}()
    H = 10 .^ range(h0, hf, length = hN)

    for h = H
        f = F21(a[1], a[2], b[1], z, h = h, order = order, taylorN = taylorN)
        push!(F, abs(f - tru))
    end
    
    return (H, F)
end

function rungridtests(path::String = ""; N = nothing)
    tests = [
        ([ .9,  1.1],   [ 1.2], 4, 99, .1, 20, 150, false, "2F1(.9, 1.1; 1.2; z)")
        ([ .9,  1.1],   [-1.2], 4, 99, .1, 20, 150, false,  "2F1(.9, 1.1; -1.2; z)")
        ([ .9, -1.1im], [-1.2], 4, 99, .1, 20, 150, false,  "2F1(.9, -1.1i; -1.2; z)")
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

function runconvergencetests(path::String = ""; N = nothing)
    tests = [
        #         a       b        z                                          tru
        ([ .9, 1.1], [ 1.2], 2 + 2im, -0.06432544451316979 + 0.5206938157025687im)
        ([-.9, 1.1], [-1.2], 2 + 2im,  2.966075098239322   + 2.865931991133645im)
        ([-.9, 1.1-2im], [-1.2im], 2 + 2im, 0.4226703551093628 - 4.373177905207582im)
    ]

    orders = [1,2,4,8,12,16,20]
    titles = ["2F1(0.9, 1.1; 1.2; 2 + 2i)", 
              "2F1(-0.9, 1.1; -1.2; 2 + 2i)",
              "2F1(-0.9, 1.1-2i; -1.2i; 2 + 2i)"]

    if isnothing(N)
        N = 1:length(tests)
    end

    h = []
    plots = Vector{Plots.Plot}()
    for (n, (name, test)) ∈ enumerate(zip(titles, tests))
        if n ∉ N
            continue
        end

        p = plot(title = name)
        for order ∈ orders
            (h, f) = convergencetest(test..., order = order)
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
        plot!(h, h, ls = :dash, lc = :black, lw = 3, label = "O(h)")
        plot!(h, h.^20, ls = :dash, lc = :gray, lw = 3, label = "O(h^20)")

        push!(plots, p)
    end
    
    if path !== ""
        for (n, p) ∈ enumerate(plots)
            savefig(p, path * "Convergence$(n).png")
        end
    end

    return plots
end
