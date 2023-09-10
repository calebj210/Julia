#=
# Collection of functions for visuallizing properties of integrators
#
# Author: Caleb Jacobs
# DLM: September 6, 2023
=#

include("Gregory.jl")
using PlotlyJS
using Colors
using BenchmarkTools

"""
    convergenceOldPlot(a, b, f, ∫f; r = 5, DFT = true, greg = false, nLow = 10, nHigh = 1000, nStep = 1, oLow = 0, oHigh = 10)
Construct convergence plot of applying DFT or Gregory corrections to `f` over the bounds from `a` to `b`. The true solution to compare to is given by `∫f`. The number of nodes away for DFT is given by `r`.
"""
function convergenceOldPlot(a, b, f, ∫f; r = 5, DFT = true, greg = false, nLow = 10, nHigh = 1000, nStep = 1, oLow = 0, oHigh = 10)
    ns   = [nLow:nStep:nHigh...]                # Number of nodes
    errsDFT     = zeros(length(ns))             # Initialize DFT errors
    errsGregory = zeros(length(ns))             # Initialize Gregory errors

    # Begin plotting
    p = plot(ns[[1,end]], 1 ./ ns[[1,end]], lw = 3, ls = :dash, legend = false)
    plot!(ns[[1,end]], 1 ./ ns[[1,end]].^2, lw = 3, ls = :dash)
    plot!(ns[[1,end]], 1 ./ ns[[1,end]].^3, lw = 3, ls = :dash)
    for N ∈ oLow : oHigh
        for i ∈ 1 : length(ns)
            if DFT
                errsDFT[i] = abs(nonSingDFTInt(f, a, b, ns[i], N, r) - ∫f)
            end

            if greg
                errsGregory[i] = abs(nonSingGregoryInt(f, a, b, ns[i], N) - ∫f)
            end
        end
        if DFT
            plot!(ns, errsDFT, 
                  yaxis = :log, 
                  xaxis = :log, 
                  lw = 2,
                  lz = N, 
                  c = :greens,
                  ylims = (1e-17, 1e0))
        end
        if greg
            plot!(ns, errsGregory, 
                  yaxis = :log, 
                  xaxis = :log, 
                  lw = 2,
                  lz = N, 
                  c = :reds,
                  ylims = (1e-17, 1e0))
        end
    end
    
    display(p)
end

"""
    convergencePlot(a, b, f, α ∫f; r = 5, DFT = true, greg = false, nLow = 10, nHigh = 1000, nStep = 1, oLow = 0, oHigh = 10)
Construct convergence plot of applying DFT or Gregory corrections to `f` over the bounds from `a` to `b`. The true solution to compare to is given by `∫f`. The number of nodes away for DFT is given by `r`.
"""
function convergencePlot(a, b, f, ∫f; α = 0, β = 0, r = 5, nLow = 10, nHigh = 1000, nStep = 1, oLow = 0, oHigh = 10)
    ns   = [nLow:nStep:nHigh...]           # Number of nodes
    errs = zeros(length(ns))               # Initialize DFT errors

    # Begin plotting
    p = plot( ns[[1,end]], 1 ./ ns[[1,end]],    lw = 3, ls = :dash, legend = false)
        plot!(ns[[1,end]], 1 ./ ns[[1,end]].^2, lw = 3, ls = :dash)
        plot!(ns[[1,end]], 1 ./ ns[[1,end]].^3, lw = 3, ls = :dash)
    for N ∈ oLow : oHigh
        for i ∈ 1 : length(ns)
            errs[i] = abs(dftInt(f, a, b, α = α, β = β, n = ns[i], N = N, r = r) - ∫f)
        end

        plot!(ns, errs, 
              yaxis = :log, 
              xaxis = :log, 
              lw = 2,
              lz = N, 
              c = :greens,
              ylims = (1e-16, 1e1))
    end
    
    display(p)
end

function complexPlot(z⃗, f⃗)
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    args = angle.(f⃗[:])
    args[args .< 0] .+= 2π

    plt = plot(heatmap(
            x = x⃗,
            y = y⃗,
            z = args,
            zsmooth = "best",
            zmin = 0, zmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["-π", "0", "π"])),
        Layout(
            width = 800, height = 800))
    
    return plt
end

function complexPlot3d(z⃗, f⃗)
    x⃗ = real(z⃗)
    y⃗ = imag(z⃗)
    z⃗ = abs.(f⃗)
    args = angle.(f⃗)
    args[args .< 0] .+= 2π

    plt = plot(surface(
            x = x⃗, y = y⃗, z = z⃗,
            surfacecolor = args,
            cmin = 0, cmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["-π", "0", "π"])))
        
    return plt
end
