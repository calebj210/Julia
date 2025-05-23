using ComplexVisuals, CairoMakie

include("2F0.jl")

function grid_test_2f0(a = 1.1, b = 1.2; tru = nothing, n = 20, m = 1)
    z = ComplexGrid(range(-50, 50, 100), range(-50, 50, 100))
    print("Expansion U... ")
    f = U.(a, b, z, n = n, m = m)
    tru_flag = false
    if isnothing(tru)
        print("done!\nMathematica U... ")
        tru = mathematica_U.(a, b, z)
        println("done!")
        
        tru_flag = true
    end

    err = abs.((f - tru) ./ tru)

    fig, ax, plt = heatmap(z, err, 
                           colorscale = log10,
                           colorrange = (1e-17, 1e0)
                          )

    Colorbar(fig[1,2], plt)
    colsize!(fig.layout, 1, Aspect(1,1))
    resize_to_layout!(fig)

    if tru_flag
        return (tru, fig)
    else
        return fig
    end
end
