using ComplexVisuals
using PlotlyJS

function loglog_plot(x::Vector{T}, y::Array{T}; targs::NamedTuple = (;), largs...) where {T <: Real}
    traces = Vector{GenericTrace}()
    for f ∈ eachcol(y)
        trace = scatter(;
            x = x, y = f,
            mode = "lines",
            targs...
        )

        push!(traces, trace)
    end

    layout = Layout(;
        plot_template_2d...,

        xaxis = attr(
            type = "log",
            autorange = "reversed"
        ),
        yaxis_type = "log",
        font_size = 15,

        largs...
    )

    plot(traces, layout)
end

function convergence_test(a, b, c, z; order = 20, taylorN = 150, h0 = -2, hf = 0, hN = 100)
    F = Vector{Float64}()
    H = 10 .^ range(h0, hf, length = hN)

    tru = mathematica_2f1(a, b, c, z)

    for h = H
        f = _2f1(a, b, c, z, H = h, order = order, N = taylorN)
        push!(F, abs(f - tru) / abs(tru))
    end
    
    return (H, F)
end
