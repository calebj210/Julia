#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: October 9, 2024
=#

using ComplexVisuals
using PlotlyJS
using MathLink
include("HypergeometricGrid.jl")

mathematica_pfq(a, b, z) = 
    try 
        Complex(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z).args...)
    catch 
        Real(weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z))
    end

"Generate pFq test"
function gridtest(a,b, r = 1, n = 20, np = 3, Tr = 0.5, modifyZ1 = true, corrR = .5, innerR = .6, outerR = .8, z1N = 70)
    if length(a) != length(b) + 1
        (z, f) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
        h = []
    else    
        (z, f, h) = pFq(a, b, r = r, n = n, np = np, Tr = Tr, modifyZ1 = modifyZ1, corrR = corrR, innerR = innerR, outerR = outerR, z1N = z1N)
    end

    tru::Vector{ComplexF64} = [mathematica_pfq(a, b, z) for z ∈ z]

    N = isqrt(length(z))
    Z = reshape(z, N, :)
    F = reshape(f, N, :)
    TRU = reshape(tru, N, :)

    title = string(length(a), "F", length(b), "(", a, "; ", b, "; z)") 

    p = generate_graphics(Z, F, TRU, title = title)

    return (z, f, h, tru, p)
end

function rungridtests(; N = 0)
    tests = 
        [
#                 a                  b                    r   n np  Tr  modz1 corR inR outR z1N
            ([1.1,1.9],          [2.9],                1.99, 41, 5, 0.5, true, .5,  .6, .8, 70),    # Test 1
            ([1.1,-1.9],         [2.9],                1.99, 41, 5, 0.5, true, .5,  .6, .8, 70),    # Test 2
            ([1.0,1.1,0.9],      [1.2,1.3],            1.99, 41, 5, 0.6, true, .5,  .6, .8, 70),    # Test 3
            ([1.0,1.1,-0.9],     [1.2,1.3],            1.99, 41, 5, 0.6, true, .5,  .6, .8, 70),    # Test 4
            ([1.0,1.1,1.2,0.9],  [1.3,1.4,1.5],        1.99, 41, 5, 0.6, true, .5,  .6, .8, 70),    # Test 5
            ([1.0,1.1,1.2,-0.9], [1.3,1.4,1.5],        1.99, 41, 5, 0.6, true, .5,  .6, .8, 70),    # Test 6
        ]

for (n, test) ∈ enumerate(tests)
        if n ∉ N
            continue
        end

        println("Running test ", test, " with a = ", a, " and b = ", b, ".")

        (z, f, h, tru, ps) = gridtest(test...)

        for p ∈ ps
            display(p)
        end

        println("Press enter when ready for next test.")
        readline()
    end

    return nothing
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

    abserr         = complex_heatmap(z, log10.(err); title_text = title)              # Absolute error plot
    relerr         = complex_heatmap(z, log10.(err ./ abs.(tru)); title_text = title) # Relative error plot
    absarg         = complex_surface_plot(z, f; title_text = title)                   # Abs-arg plot
    replot, implot = complex_reim_surface_plot(z, f, title_text = title)              # Real and imaginary part plots

    return (abserr, relerr, absarg, replot, implot)
end

generate_graphics(z::T, f::T, tru::T; title = "") where T <: Vector{<: Number} = 
begin
    N = isqrt(length(z))

    Z   = reshape(z,   N, N)
    F   = reshape(f,   N, N)
    TRU = reshape(tru, N, N)

    generate_graphics(Z, F, TRU; title = title)
end
    
