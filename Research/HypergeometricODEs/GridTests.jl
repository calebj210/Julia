using ComplexVisuals
using PlotlyJS
include("pFq.jl")

function gridtest(a, b, c, r, n, H = .1, order = 20, taylorN = 100, title = "")
    z = complex_grid(r, n);

    tru = mathematica_2f1.(a, b, c, z)
    val = _2f1.(a, b, c, z)
    
    plts = generate_graphics(z, val, tru, title = title)
    
    return plts
end

function rungridtests(path::String = ""; N = nothing)
    tests = [
        ([-.9,  1.11], [1.2],    4, 200, .1,  150, 150, false, "2F1(-.9, 1.11; 1.2; z)")

        ([-.9, 5.0-20im], [1.2], 4, 200, .05, 150, 150, false, "2F1(-.9, 5-20im; 1.2; z)")
        ([-.9, 1.11], [50.0im],  20, 200, .05, 150, 150, false, "2F1(-.9, 1.11; 50i; z)")
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

"Generate a heatmap over a complex domain"
function complex_heatmap(z::T, f::S; kwargs...) where {T <: Vector{<: Number}, S <: Vector{<: Real}}
    x, y = reim(z)
    
    trace = heatmap(
        x = x, y = y, z = f,

        colorscale = "Viridis",
        zmin = -17, zmax = 1,
    )

    layout = Layout(;
        plot_template_2d_complex...,
        kwargs...
    )
    
    return plot(trace, layout)
end
complex_heatmap(z::T, f::S; kwargs...) where {T <: Matrix{<: Number}, S <: Matrix{<: Real}} = complex_heatmap(vec(z), vec(f); kwargs...)

"Generate graphics"
function generate_graphics(z::T, f::T, tru::T; title = "") where T <: Matrix{<: Number}
    err = abs.(f - tru)
    err[iszero.(err)] .= 1e-17

    abserr         = complex_heatmap(z, log10.(err))                # Absolute error plot
    relerr         = complex_heatmap(z, log10.(err ./ abs.(tru)))   # Relative error plot
    absarg         = complex_surface_plot(z, f; title_text = title) # Abs-arg plot
    replot, implot = complex_reim_surface_plot(z, f)                # Real and imaginary part plots

    return (abserr, relerr, absarg, replot, implot)
end
