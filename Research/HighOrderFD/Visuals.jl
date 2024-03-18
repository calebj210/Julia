#=
# Plotting routines for FD stencils
#
# Author: Caleb Jacobs
# DLM: March 17, 2024
=#
using Plots
using LaTeXStrings

"""
    distancePlot2D(Z, D)
Create a log plot of the absolute value of the stencil weights in `D` as a function of distance from the origin for nodes in `Z`.
"""
function distancePlot2D(Z, D; title = "", n = 2)
    ω = vec(Float64.(abs.(D)))
    z = vec(norm.(Z))

    ω[ω .== 0] .= ω[1]

    x = range(0, maximum(z), length = 100)
    y = maximum(ω) * exp.(-π / n * x.^2)

    p = plot(x, y,
         ls = :dash,
         lc = :black,

         yscale = :log10,
         ylims  = (ω[1], 20maximum(ω)),
         title  = title,
         xlabel = "Distance from Origin",
         ylabel = "Abs-Weight",
         label  = latexstring(raw"C e^{-\frac{\pi}{" * string(n) * raw"}x^2}"),
         
         dpi = 300)


    plot!(z, ω, 
         st  = :scatter,
         mc  = :gray,
         mα  = .4,
         msw = 0,

         label = "Weights")

    return p
end
