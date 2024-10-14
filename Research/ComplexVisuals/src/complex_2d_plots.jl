using PlotlyJS
using Images, Base64
const CS = PlotlyJS.PlotlyBase.ColorSchemes

const plot_template_2d = (;
    template = "plotly_white",
    width = 1000, height = 1000,

    title = attr(
        font_size = 25,
        y = 0.96,
        x = 0.5
    ),

    font_size = 20,

    xaxis_title_font_size = 25,
    yaxis_title_font_size = 25,
)

const plot_template_2d_complex = (;
    plot_template_2d...,

    xaxis_title_text = "Re(z)",
    yaxis_title_text = "Im(z)",
)

function complex_phase_plot(z::T, f::T; kwargs...) where T <: Vector{<: Number}
    x, y = reim(z)
    phase = mod.(angle.(f), 2π)

    trace = heatmap(
        x = x,
        y = y,
        z = phase,

        zsmooth = "fast",
        colorscale = CS.hsv,
        zmin = 0, zmax = 2π,
        showscale = false,
    )

    png_data = read("$(@__DIR__)/AngleWheel.png", String) # Read the PNG file as a string
    png_base64 = base64encode(png_data)

    layout = Layout(;
        plot_template_2d_complex...,

        images = [attr(
            source = "data:image/png;base64,$png_base64",
            xref = "paper",  # Reference the paper coordinates (0 to 1)
            yref = "paper",
            x = 0.8,         # Position horizontally (0.5 is center)
            y = 0.2,         # Position vertically
            sizex = 0.2,      # Size relative to paper width
            sizey = 0.2,      # Size relative to paper height
            opacity = 0.8,    # Adjust transparency as needed
            layer = "above"   # Place the image above the 3D graph
        )],

        kwargs...
    )

    plot(trace, layout)
end
complex_phase_plot(z::T, f::T; kwargs...) where T <: Matrix{<: Number} = complex_phase_plot(vec(z), vec(f); kwargs...)
