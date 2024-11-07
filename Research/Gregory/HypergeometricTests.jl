#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: November 5, 2024
=#

# using ComplexVisuals
using ComplexVisualsMakie, CairoMakie
# using PlotlyJS
using MathLink
using HypergeometricFunctions: pFqweniger, pFqdrummond
using Nemo: hypergeometric_1f1, hypergeometric_2f1, ComplexField # Arb/Flint interface for Julia
using Suppressor
include("HypergeometricGrid.jl")

function jo_2f1(a, b, c, z::Number) :: ComplexF64
    cc = ComplexField()
    (a,b,c,z) = cc.((a,b,c,z))

    return hypergeometric_2f1(a, b, c, z)
end

mathematica_pfq(a, b, z) = 
    try 
        Complex(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z).args...)
    catch 
        Real(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z))
    end

"Generate pFq test"
function gridtest(a,b; kwargs...)
        defaults = (;
            grid_radius = 1, 
            grid_points = 41, 
            padding_layers = 5, 
            taylor_radius = 0.5, 
            correction_radius = .5, 
            inner_radius = .6, 
            outer_radius = .8, 
            z1_expansion_order = 70,
            modify_z1 = true) 
        kwargs = merge(defaults, kwargs)

    if length(a) != length(b) + 1
        (z, f) = pFq(a, b; kwargs...)
        h = []
    else    
        println("Tra time: ")
        @time (z, f, h) = pFq(a, b; kwargs...)
    end

    println("Wen time: ")
    @time wenf = @suppress_err [pFqweniger(a, b, z) for z ∈ z]

    # println("Dru time: ")
    # @time druf = @suppress_err [pFqdrummond(a, b, z) for z ∈ z]

    println("Jo time: ")
    @time jof = @suppress_err jo_2f1.(a..., b..., z)

    println("Mat time: ")
    @time tru::Vector{ComplexF64} = [mathematica_pfq(a, b, z) for z ∈ z]

    title = string(length(a), "F", length(b), "(", a, "; ", b, "; z)") 
    
    trap_plots = generate_graphics(z, f, tru, title = title)
    wen_plots = generate_graphics(z, wenf, tru, title = title)
    # dru_plots = generate_graphics(z, druf, tru, title = title)
    jo_plots = generate_graphics(z, jof, tru, title = title)
    
    return (;trap = trap_plots, wen = wen_plots, jo = jo_plots) #, dru = dru_plots)
end

function large_b_test()
    a = [1.1, 6.32]
    b = [1.2]
    
    test = (;
            grid_radius = 2.49,
            grid_points = 71,
            padding_layers = 5,
            taylor_radius = .5,
            correction_radius = .7,
            inner_radius = .7,
            outer_radius = .8,
            z1_expansion_order = 120,
            modify_z1 = true
           )
    
    results = gridtest(a, b; test...)
end

function pattern_test()
    a = [1.99, 0.9]
    b = [2.9]
    
    test = (;
            grid_radius = 2.49,
            grid_points = 81,
            padding_layers = 5,
            taylor_radius = .5,
            correction_radius = .25,
            inner_radius = .6,
            outer_radius = .8,
            z1_expansion_order = 120,
            modify_z1 = true
           )

    (z,f,h) = pFq(a, b; test...)

    tru = [mathematica_pfq(a, b, z) for z ∈ z]

    relerr = abs.((f - tru) ./ tru)

    fig,ax,plt = heatmap(real(z), imag(z), relerr;
                         colorscale = log10,
                         colorrange = (1e-17, 1e-6))
    lines!(ax, [-2.49, 2.49], [0, 0], linewidth = 1, color = :white)
    lines!(ax, [0, 0], [-2.49, 2.49], linewidth = 1, color = :white)
    lines!(ax, reim(cispi.(range(0, 2π, 100)))..., linewidth = 1, color = :white)
    ax.xlabel = "Re(z)"
    ax.ylabel = "Im(z)"
    ax.title  = L"End-Corrected Trapezoidal Rule for ${_2}F_1(1.99, 0.9; 2.9; z)$"
    Colorbar(fig[1,2], plt)
    colsize!(fig.layout, 1, Aspect(1, 1))
    resize_to_layout!(fig)

    return Makie.FigureAxisPlot(fig, ax, plt)
end

function grid_time_complexity_error(a, b, c; r = 2, Ns = 81:22:301)
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

    fig, ax, plt = lines(Ns, times)

    return fig
end

"Generate a heatmap over a complex domain"
# function complex_heatmap(z::T, f::S; kwargs...) where {T <: Vector{<: Number}, S <: Vector{<: Real}}
#     x, y = reim(z)
    
#     trace = heatmap(
#         x = x, y = y, z = f,

#         colorscale = "Viridis",
#         zmin = -17, zmax = 1,
#     )

#     layout = Layout(;
#         plot_template_2d_complex...,
#         kwargs...
#     )
    
#     return plot(trace, layout)
# end
# complex_heatmap(z::T, f::S; kwargs...) where {T <: Matrix{<: Number}, S <: Matrix{<: Real}} = complex_heatmap(vec(z), vec(f); kwargs...)

"Generate graphics"
# function generate_graphics(z::T, f::T, tru::T; title = "") where T <: Matrix{<: Number}
#     err = abs.(f - tru)
#     err[iszero.(err)] .= 1e-17

#     abserr         = complex_heatmap(z, log10.(err); title_text = title)              # Absolute error plot
#     relerr         = complex_heatmap(z, log10.(err ./ abs.(tru)); title_text = title) # Relative error plot

#     f[imag(z) .≈ 0 .&& real(z) .>= 1] .= NaN + NaN * im

#     absarg         = complex_surface_plot(z, f; title_text = "Abs-Arg(f)")            # Abs-arg plot
#     replot, implot = complex_reim_surface_plot(z, f)                                  # Real and imaginary part plots

#     return (abserr, relerr, absarg, replot, implot)
# end

# generate_graphics(z::T, f::T, tru::T; title = "") where T <: Vector{<: Number} = 
# begin
#     N = isqrt(length(z))

#     Z   = reshape(z,   N, N)
#     F   = reshape(f,   N, N)
#     TRU = reshape(tru, N, N)

#     generate_graphics(Z, F, TRU; title = title)
# end
